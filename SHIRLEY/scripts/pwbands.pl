#!/usr/bin/perl -w

sub klen {
  ($kx, $ky, $kz)=@_;
  sqrt($kx*$kx + $ky*$ky + $kz*$kz);
}

sub kdiff {
  ($kx1, $ky1, $kz1, $kx2, $ky2, $kz2) = @_;
  ($kx2-$kx1, $ky2-$ky1, $kz2-$kz1);
}

# find alat first
while(<STDIN>) {
  chomp;
  if( m/lattice parameter/ ) {
    @data1=split /=/;
    @data2=split ' ', $data1[1];
    $alat=$data2[0];
    last;
  }
}

$prefac=2.0*3.141592654/$alat;
print "# lattice parameter = $alat\n";
print "#         prefactor = $prefac\n";

$spin=0;
while(<STDIN>) {
#  if( m/End of band structure calculation/ ) {
  if( m/End of .* calculation/ ) {
    $firstk=1;
    <STDIN>;
    $_=<STDIN>;
    while( m/SPIN/ ) {
      $spin++;
      print "#spin $spin\n";
      <STDIN>;
      <STDIN>;
      chomp($line=<STDIN>);
      $kpt_count=0;
      while( $line=~/k =/ ) {
        @data1 = split /k =/, $line;
        @data2 = split /  /, $data1[1];
        @kpt = ( substr($data2[0], 0,7), 
                 substr($data2[0], 7,7), 
                 substr($data2[0],14,7) );
        $kpt_count++;
        $eig='';
        $_=<STDIN>;
        chomp($line=<STDIN>);
        until($line =~ /^ *$/) {
          $eig.=$line;
          chomp($line=<STDIN>);
        }
        if( $firstk ) {
          $kpt_pathlen=0.0;
          undef($firstk);
        } else {
          @dkpt = &kdiff( @kpt_last, @kpt );
          $kpt_pathlen += &klen(@dkpt)*$prefac;
        }
        @kpt_last = @kpt;
        print "$kpt_count @kpt $kpt_pathlen $eig\n";
        chomp($line=<STDIN>);
      }
      $_=$line;
      $firstk=1;
    } 
    if($spin>0) { exit; }

    chomp($line=$_);
    $kpt_count=0;
    while( $line=~/k =/ ) {
      @data1 = split /k =/, $line;
      @data2 = split /  /, $data1[1];
      @kpt = ( substr($data2[0], 0,7), 
               substr($data2[0], 7,7), 
               substr($data2[0],14,7) );
      $kpt_count++;
      $eig='';
      <STDIN>;
      chomp($line=<STDIN>);
      until($line =~ /^ *$/) {
#        print "line: $line\n";
        $eig.=$line;
        chomp($line=<STDIN>);
#        print "new line: $line\n";
      }
      if( $firstk ) {
        $kpt_pathlen=0.0;
        undef($firstk);
      } else {
        @dkpt = &kdiff( @kpt_last, @kpt );
        $kpt_pathlen += &klen(@dkpt)*$prefac;
      }
      @kpt_last = @kpt;
      print "$kpt_count @kpt $kpt_pathlen $eig\n";
      chomp($line=<STDIN>);
    }
  }
}
exit;

