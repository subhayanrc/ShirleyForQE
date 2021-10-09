#!/usr/bin/perl -w

use Getopt::Long;

$usage=" pwstress.pl
  Syntax: pwstress.pl [option]
  Reads input from STDIN and looks for stress info which is dumped 
  as a single line depending on the options provided:
    -help|-h  print this message
    -P        print total pressure (kbar)
    -p        print pressure tensor (kbar)
    -s        print stress tensor (ryd/bohr**3)
  Ordering of stress tensors is
    1=xx ; 2=yy ; 3=zz ; 4=yz ; 5=xz ; 6=xy
";

Getopt::Long::Configure ("bundling");
GetOptions(
   'help|h'      => \$help,
   'P'           => \$show_presstot,
   'p'           => \$show_press,
   's'           => \$show_stress,
);

if( $help ) {
  print $usage;
  exit;
}

# Imagine some output like this:
#    entering subroutine stress ...
#
#       total   stress  (ryd/bohr**3)                  (kbar)     P=    0.00
#     0.00000003   0.00000000   0.00000000          0.00      0.00      0.00
#     0.00000000   0.00000003   0.00000000          0.00      0.00      0.00
#     0.00000000   0.00000000  -0.00000002          0.00      0.00      0.00
#
while($_=<STDIN>) {
  if( m/entering subroutine stress .../ ) {
    chomp($_=<STDIN>);

    chomp($_=<STDIN>);
    @data=split /P=/, $_;
    $presstot=$data[1];

    # xx, xy, xz
    chomp($_=<STDIN>);
    @data=split ' ', $_;
    $stress[0]=$data[0];
    $stress[5]=$data[1];
    $stress[4]=$data[2];
    $press[0]=$data[3];
    $press[5]=$data[4];
    $press[4]=$data[5];

    # yy, yz
    chomp($_=<STDIN>);
    @data=split ' ', $_;
    $stress[1]=$data[1];
    $stress[3]=$data[2];
    $press[1]=$data[4];
    $press[3]=$data[5];

    # zz
    chomp($_=<STDIN>);
    @data=split ' ', $_;
    $stress[2]=$data[2];
    $press[2]=$data[5];

    if( $show_presstot ) {
      printf "%f\n", $presstot;
    } elsif ( $show_press ) {
      printf "%f %f %f %f %f %f\n", @press;
    } elsif ( $show_stress ) {
      printf "%f %f %f %f %f %f\n", @stress;
    } else {
      printf "%f %f %f %f %f %f %f %f %f %f %f %f %f\n", $presstot, @press, @stress;
    }

  }
}
