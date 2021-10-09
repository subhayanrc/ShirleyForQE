#!/usr/bin/perl -w

# find nelec first
while(<STDIN>) {
  chomp;
  if( m/number of electrons *= *(.*) *$/ ) {
    $nelec=$1;
    last;
  }
}

print "# nelec = $nelec\n";
$homo=int($nelec/2);
print "# HOMO is state $homo\n";
$homo--;

while(<STDIN>) {
  chomp;
  if( m/End of band structure calculation/ ||
      m/ End of self-consistent calculation/ ) {
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
      <STDIN>;
      chomp($line=<STDIN>);
      until($line =~ /^ *$/) {
#        print "line: $line\n";
        $eig.=$line;
        chomp($line=<STDIN>);
#        print "new line: $line\n";
      }
      print "$kpt_count @kpt $eig\n";
      @eigval=split ' ', $eig;
      if( $kpt_count > 1 ) {
        print "# More than one k-point - HOMO may not be meaningful\n";
      }
      print "# HOMO energy = $eigval[$homo]\n";
      chomp($line=<STDIN>);
    }
  }
}
exit;

