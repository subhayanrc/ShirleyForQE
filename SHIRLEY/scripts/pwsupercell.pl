#!/usr/bin/perl

use POSIX;

# pwsupercell.pl
$usage = " usage: pwsupercell.pl N M L
                  to create a supercell of NxMxL unit cells
                  with input read from STDIN
";

die "$usage" if @ARGV != 3;
$N = shift @ARGV;
$M = shift @ARGV;
$L = shift @ARGV;

# read the input file
@infile=();
while(<STDIN>) {
  push @infile, $_;
}

# matching regex for a number (integer, decimal, or float)
$num=qr/[+-]?(\d*\.?\d+|\d+\.?\d*)([eE][+-]?\d+)?/;
# find details
foreach (@infile) {
  if(/\bnat *= *(\d+)/) {
    $nat=$1;
  } 
  if(/\bnbnd *= *(\d+)/) {
    $nbnd=$1;
  } 
  if(/\bibrav *= *(\d+)/) {
    $ibrav=$1;
  } 
  if(/\bcelldm\(1\) *= *($num)/) {
    $celldm[1]=$1;
  } 
  if(/\bcelldm\(2\) *= *($num)/) {
    $celldm[2]=$1;
  } 
  if(/\bcelldm\(3\) *= *($num)/) {
    $celldm[3]=$1;
  } 
  if(/\bcelldm\(4\) *= *($num)/) {
    $celldm[4]=$1;
  } 
  if(/\bcelldm\(5\) *= *($num)/) {
    $celldm[5]=$1;
  } 
  if(/\bcelldm\(6\) *= *($num)/) {
    $celldm[6]=$1;
  } 
  if(/\ba *= *($num)/) {
    $a=$1;
  } 
  if(/\bb *= *($num)/) {
    $b=$1;
  } 
  if(/\bc *= *($num)/) {
    $c=$1;
  } 
  if(/\bcosab *= *($num)/) {
    $cosab=$1;
  } 
  if(/\bcosac *= *($num)/) {
    $cosac=$1;
  } 
  if(/\bcosbc *= *($num)/) {
    $cosbc=$1;
  } 
}

#print "nat=$nat\n";
#print "ibrav=$ibrav\n";
#if($celldm[1]) {
#  print "celldm(1) = $celldm[1]\n";
#  if($celldm[2]) { print "celldm(2) = $celldm[2]\n" };
#  if($celldm[3]) { print "celldm(3) = $celldm[3]\n" };
#  if($celldm[4]) { print "celldm(4) = $celldm[4]\n" };
#  if($celldm[5]) { print "celldm(5) = $celldm[5]\n" };
#  if($celldm[6]) { print "celldm(6) = $celldm[6]\n" };
#} else {
#  print "a=$a\n";
#  if($b) { print "b=$b\n" };
#  if($c) { print "c=$c\n" };
#  if($cosab) { print "cosab=$cosab\n" };
#  if($cosac) { print "cosac=$cosac\n" };
#  if($cosbc) { print "cosbc=$cosbc\n" };
#}

# find lattice
$n=-1;
foreach(@infile) {
  $n++;
  if(/CELL_PARAMETERS *\{* *(\w*)/) {
    if($1) { $cubhex=$1 }
    $param_offset=$n;
    for($i=1; $i<=3; $i++) {
      @data=split ' ', $infile[$n+$i];
      $at[1][$i]=$data[0];
      $at[2][$i]=$data[1];
      $at[3][$i]=$data[2];
      #print "$data[0] $data[1] $data[2] -> $at[1][$i] $at[2][$i] $at[3][$i]\n";
    }
    # remove cell params now that we have them
    splice(@infile,$param_offset,4);
  }
}

# find atoms
$n=-1;
foreach (@infile) {
  $n++;
  if(/ATOMIC_POSITIONS *[\{(]* *(\w+)/) {
    $unit=$1;
    for($i=1; $i<=$nat; $i++) {
      @data=split ' ', $infile[$n+$i];
      $label[$i]=$data[0];
      $atom[1][$i]=$data[1];
      $atom[2][$i]=$data[2];
      $atom[3][$i]=$data[3];
    }
  }
}
#print "unit = $unit\n";
#for($i=1; $i<=$nat; $i++) {
#  print "$label[$i] $atom[1][$i] $atom[2][$i] $atom[3][$i]\n";
#}

# make lattice matrix
$bohrang=0.529177;
if($a) { $celldm[1]=$a/$bohrang; }
if($b) { $celldm[2]=$b/$a };
if($c) { $celldm[3]=$c/$a };
if($ibrav!=14) {
  if($cosab) { $celldm[4]=$cosab }
} else {
  if($cosbc) { $celldm[4]=$cosbc } # alpha
  if($cosac) { $celldm[5]=$cosac } # beta
  if($cosab) { $celldm[6]=$cosab } # gamma
}
if($ibrav==0 && $celldm[1]) {
  $alat=$celldm[1];
} elsif($ibrav==0 && ! $celldm[1]) {
  $celldm[1]=sqrt( $at[1][1]*$at[1][1] + $at[1][2]*$at[1][2] + $at[1][3]*$at[1][3] );
  $alat=$celldm[1];
  for($j=1; $j<=3; $j++) {
    for($i=1; $i<=3; $i++) {
      $at[$i][$j] /= $alat;
    }
  }
} else {
# generate lattice
  $at[1][1]=0.0;
  $at[1][2]=0.0;
  $at[1][3]=0.0;
  $at[2][1]=0.0;
  $at[2][2]=0.0;
  $at[2][3]=0.0;
  $at[3][1]=0.0;
  $at[3][2]=0.0;
  $at[3][3]=0.0;
  if ($ibrav==1) {
    # simple cubic
    $at[1][1]=$celldm[1];
    $at[2][2]=$celldm[1];
    $at[3][3]=$celldm[1];
  } elsif($ibrav==2) {
    # fcc
    $term=$celldm[1]*0.5;
    $at[1][1]=-$term;
    $at[1][3]= $term;
    $at[2][2]= $term;
    $at[2][3]= $term;
    $at[3][1]=-$term;
    $at[3][2]= $term;
  } elsif($ibrav==3) {
    # bcc
    $term=$celldm[1]*0.5;
    $at[1][1]= $term;
    $at[1][2]= $term;
    $at[1][3]= $term;
    $at[2][1]=-$term;
    $at[2][2]= $term;
    $at[2][3]= $term;
    $at[3][1]=-$term;
    $at[3][2]=-$term;
    $at[3][3]= $term;
  } elsif($ibrav==4) {
    # hexagonal
    $cbya=$celldm[3];
    $at[1][1]=$celldm[1];
    $at[2][1]=-$celldm[1]*0.5;
    $at[2][2]=$celldm[1]*sqrt(3.0)*0.5;
    $at[3][3]=$celldm[1]*$cbya;
  } elsif($ibrav==5) {
    # trigonal
    $term1=sqrt(1.0+2.0*$celldm[4]);
    $term2=sqrt(1.0-$celldm[4]);
    $at[2][2]=sqrt(2.0)*$celldm[1]*$term2/sqrt(3.0);
    $at[2][3]=$celldm[1]*$term1/sqrt(3.0);
    $at[1][1]=$celldm[1]*$term2/sqrt(2.0);
    $at[1][2]=-$at[1][1]/sqrt(3.0);
    $at[1][3]= $at[2][3];
    $at[3][1]=-$at[1][1];
    $at[3][2]= $at[1][2];
    $at[3][3]= $at[2][3];
  } elsif($ibrav==6) {
    # tetragonal
    $cbya=$celldm[3];
    $at[1][1]=$celldm[1];
    $at[2][2]=$celldm[1];
    $at[3][3]=$celldm[1]*$cbya;
  } elsif($ibrav==7) {
    # body-centered tetragonal
    $cbya=$celldm[3];
    $at[2][1]=$celldm[1]*0.5;
    $at[2][2]=$at[2][1];
    $at[2][3]=$cbya*$celldm[1]*0.5;
    $at[1][1]= $at[2][1];
    $at[1][2]=-$at[2][1];
    $at[1][3]= $at[2][3];
    $at[3][1]=-$at[2][1];
    $at[3][2]=-$at[2][1];
    $at[3][3]= $at[2][3];
  } elsif($ibrav==8) {
    # orthorhombic
    $at[1][1]=$celldm[1];
    $at[2][2]=$celldm[1]*$celldm[2];
    $at[3][3]=$celldm[1]*$celldm[3];
  } elsif($ibrav==9) {
    # one face-centered orthorhombic
    $at[1][1]= 0.5*$celldm[1];
    $at[1][2]= $at[1][1]*$celldm[2];
    $at[2][1]=-$at[1][1];
    $at[2][2]= $at[1][2];
    $at[3][3]= $celldm[1]*$celldm[3];
  } elsif($ibrav==10) {
    # all face-centered orthorhombic
    $at[2][1]= 0.5*$celldm[1];
    $at[2][2]= $at[2][1]*$celldm[2];
    $at[1][1]= $at[2][1];
    $at[1][3]= $at[2][1]*$celldm[3];
    $at[3][2]= $at[2][1]*$celldm[2];
    $at[3][3]= $at[1][3];
  } elsif($ibrav==11) {
    # body-centered orthorhombic
    $at[1][1]= 0.5*$celldm[1];
    $at[1][2]= $at[1][1]*$celldm[2];
    $at[1][3]= $at[1][1]*$celldm[3];
    $at[2][1]=-$at[1][1];
    $at[2][2]= $at[1][2];
    $at[2][3]= $at[1][3];
    $at[3][1]=-$at[1][1];
    $at[3][2]=-$at[1][2];
    $at[3][3]= $at[1][3];
  } elsif($ibrav==12) {
    # monoclinic
    $sen=sqrt(1.0-$celldm[4]*$celldm[4]);
    $at[1][1]=$celldm[1];
    $at[2][1]=$celldm[1]*$celldm[2]*$celldm[4];
    $at[2][2]=$celldm[1]*$celldm[2]*$sen;
    $at[3][3]=$celldm[1]*$celldm[3];
  } elsif($ibrav==13) {
    # one face-centered monoclinic
    $sen=sqrt(1.0-$celldm[4]*$celldm[4]);
    $at[1][1]=0.5*$celldm[1]*$sen;
    $at[1][2]=0.5*$celldm[1]*($celldm[4]-$celldm[2]);
    $at[2][1]=0.5*$celldm[1]*$sen;
    $at[2][2]=0.5*$celldm[1]*($celldm[4]+$celldm[2]);
    $at[3][3]=$celldm[1]*$celldm[3];
  } elsif($ibrav==14) {
    # triclinic
    $singam=sqrt(1.0-$celldm[6]*$celldm[6]);
    $term  =sqrt((1.0+2.0*$celldm[4]*$celldm[5]*$celldm[6]
                         -$celldm[4]*$celldm[4]
                         -$celldm[5]*$celldm[5]
                         -$celldm[6]*$celldm[6])/(1.0-$celldm[6]*$celldm[6]));
    $at[1][1]=$celldm[1];
    $at[2][1]=$celldm[1]*$celldm[2]*$celldm[6];
    $at[2][2]=$celldm[1]*$celldm[2]*$singam;
    $at[3][1]=$celldm[1]*$celldm[3]*$celldm[5];
    $at[3][2]=$celldm[1]*$celldm[3]*($celldm[4]-$celldm[5]*$celldm[6])/$singam;
    $at[3][3]=$celldm[1]*$celldm[3]*$term;
  } else {
    die "nonexistent bravais lattice: $ibrav\n";
  }
  $alat=$celldm[1];
  $at[1][1]/=$alat;
  $at[1][2]/=$alat;
  $at[1][3]/=$alat;
  $at[2][1]/=$alat;
  $at[2][2]/=$alat;
  $at[2][3]/=$alat;
  $at[3][1]/=$alat;
  $at[3][2]/=$alat;
  $at[3][3]/=$alat;

  $a=$alat*$bohrang;
}
print " alat = $alat\n";
print " at = \n";
print " $at[1][1] $at[2][1] $at[3][1]\n";
print " $at[1][2] $at[2][2] $at[3][2]\n";
print " $at[1][3] $at[2][3] $at[3][3]\n";

# convert coordinates to crystal lattice
if($unit =~ /angstrom/) {
  # make the inverse transformation
  $det=$at[1][1]*($at[3][3]*$at[2][2]-$at[3][2]*$at[2][3])
      -$at[2][1]*($at[3][3]*$at[1][2]-$at[3][2]*$at[1][3])
      +$at[3][1]*($at[2][3]*$at[1][2]-$at[2][2]*$at[1][3]);
  $det*=$alat*$bohrang;
  $bg[1][1]= ($at[3][3]*$at[2][2]-$at[3][2]*$at[2][3])/$det;
  $bg[2][1]=-($at[3][3]*$at[2][1]-$at[3][1]*$at[2][3])/$det;
  $bg[3][1]= ($at[3][2]*$at[2][1]-$at[3][1]*$at[2][2])/$det;
  $bg[1][2]=-($at[3][3]*$at[1][2]-$at[3][2]*$at[1][3])/$det;
  $bg[2][2]= ($at[3][3]*$at[1][1]-$at[3][1]*$at[1][3])/$det;
  $bg[3][2]=-($at[3][2]*$at[1][1]-$at[3][1]*$at[1][2])/$det;
  $bg[1][3]= ($at[2][3]*$at[1][2]-$at[2][2]*$at[1][3])/$det;
  $bg[2][3]=-($at[2][3]*$at[1][1]-$at[2][1]*$at[1][3])/$det;
  $bg[3][3]= ($at[2][2]*$at[1][1]-$at[2][1]*$at[1][2])/$det;

#print " bg = \n";
#print " $bg[1][1] $bg[2][1] $bg[3][1]\n";
#print " $bg[1][2] $bg[2][2] $bg[3][2]\n";
#print " $bg[1][3] $bg[2][3] $bg[3][3]\n";

  # convert coordinates to crystal
  for ($i=1; $i<=$nat; $i++) {
    for ($j=1; $j<=3; $j++) {
      $catom[$j] = $bg[$j][1]*$atom[1][$i]
                 + $bg[$j][2]*$atom[2][$i]
                 + $bg[$j][3]*$atom[3][$i];
    }
    $atom[1][$i]=$catom[1];
    $atom[2][$i]=$catom[2];
    $atom[3][$i]=$catom[3];
    print "$label[$i] $atom[1][$i] $atom[2][$i] $atom[3][$i]\n";
  }
}

# finally make the supercell
$pos="";
for($n=0; $n<$N; $n++) {
for($m=0; $m<$M; $m++) {
for($l=0; $l<$L; $l++) {
  $s[1]=(1.0/$N);
  $s[2]=(1.0/$M);
  $s[3]=(1.0/$L);

  $t[1]=$n*(1.0/$N);
  $t[2]=$m*(1.0/$M);
  $t[3]=$l*(1.0/$L);

  for($i=1; $i<=$nat; $i++) {
    $catom[1]=$atom[1][$i]*$s[1]+$t[1];
    $catom[2]=$atom[2][$i]*$s[2]+$t[2];
    $catom[3]=$atom[3][$i]*$s[3]+$t[3];

    if($unit =~ /angstrom/) {
      # convert to Cartesian
      for($j=1; $j<=3; $j++) {
        $datom[$j]=$N*$at[$j][1]*$catom[1]
                  +$M*$at[$j][2]*$catom[2]
                  +$L*$at[$j][3]*$catom[3];
        $datom[$j]*=$alat*$bohrang;
      }
      $pos=$pos . "$label[$i] $datom[1] $datom[2] $datom[3]\n";
    } else {
      $pos=$pos . "$label[$i] $catom[1] $catom[2] $catom[3]\n";
    }
  }
}
}
}

# remove old atomic positions and replace them with new
$n=-1;
foreach (@infile) {
  $n++;
  if(/ATOMIC_POSITIONS/) { $atom_offset=$n }
}
splice(@infile,$atom_offset+1,$nat,$pos);

#update k-points
$n=-1;
foreach (@infile) {
  $n++;
  if(/K_POINTS.*automatic/) { 
    @data=split ' ', $infile[$n+1];
    $nk1=shift @data;
    $nk2=shift @data;
    $nk3=shift @data;

    $nk1=ceil($nk1/$N); if($nk1<1) { $nk1=1 } 
    $nk2=ceil($nk2/$M); if($nk2<1) { $nk2=1 }
    $nk3=ceil($nk3/$L); if($nk3<1) { $nk3=1 }

    #update
    $infile[$n+1]=" $nk1 $nk2 $nk3 @data\n";
  }
}


# celldm[1]
#$celldm[1]=sqrt($at[1][1]*$at[1][1]
#               +$at[1][2]*$at[1][2]
#               +$at[1][3]*$at[1][3]);
# rescale parameters
for($i=1; $i<=3; $i++) {
  $at[$i][1]=$at[$i][1]*$N;
  $at[$i][2]=$at[$i][2]*$M;
  $at[$i][3]=$at[$i][3]*$L;
}
#$celldm[1]=sqrt($at[1][1]*$at[1][1]
#               +$at[1][2]*$at[1][2]
#               +$at[1][3]*$at[1][3]) /
#           $celldm[1];
#$celldm[1]*=$alat;

# update cell dimensions
$nat=$nat*$N*$M*$L;
$nbnd=$nbnd*$N*$M*$L;
$n=-1;
foreach (@infile) {
  $n++;
  $infile[$n] =~ s/\bnat *= *$num/nat=$nat/;
  $infile[$n] =~ s/\bnbnd *= *$num/nbnd=$nbnd/;
  $infile[$n] =~ s/\bibrav *= *$num/ibrav=0/;
  $infile[$n] =~ s/\ba *= *$num/a=$a/;
  $infile[$n] =~ s/\bb *= *$num//;
  $infile[$n] =~ s/\bc *= *$num//;
  $infile[$n] =~ s/\bcosbc *= *$num//;
  $infile[$n] =~ s/\bcosac *= *$num//;
  $infile[$n] =~ s/\bcosab *= *$num//;
  $infile[$n] =~ s/\bcelldm\(1\) *= *$num/celldm\(1\)=$alat/;
  $infile[$n] =~ s/\bcelldm\(2\) *= *$num//;
  $infile[$n] =~ s/\bcelldm\(3\) *= *$num//;
  $infile[$n] =~ s/\bcelldm\(4\) *= *$num//;
  $infile[$n] =~ s/\bcelldm\(5\) *= *$num//;
  $infile[$n] =~ s/\bcelldm\(6\) *= *$num//;
}

$param=sprintf(" CELL_PARAMETERS
  %f %f %f
  %f %f %f
  %f %f %f\n",
  $at[1][1],$at[2][1],$at[3][1],
  $at[1][2],$at[2][2],$at[3][2],
  $at[1][3],$at[2][3],$at[3][3]);
push @infile, $param;

# now dump the modified input file
print @infile;
exit;
