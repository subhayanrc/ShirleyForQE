#!/usr/bin/perl -w

use Getopt::Long;

$usage=" pwforce.pl
  Syntax: pwforce.pl [option]
  Reads input from STDIN and looks for stress info which is dumped 
  as a single line depending on the options provided:
    -help|-h        print this message
    -F              print total force
    [-a,--atom=]a   print force on atom a
    [-t,--type=]t   print forces on atom of type t
  Format of output is
    Type F_x F_y F_z
";

Getopt::Long::Configure ("bundling");
GetOptions(
   'help|h'      => \$help,
   'a|atom=i'      => \$show_atom,
   't|type=i'      => \$show_type,
   'F'           => \$show_forcetot
);

if( $help ) {
  print $usage;
  exit;
}

#     Forces acting on atoms (Ry/au):
#
#     atom   1 type  1   force =     0.00000000    0.00000000    0.00000000
#     atom   2 type  2   force =     0.00000000   -0.00000015    0.00000000
#     atom   3 type  2   force =     0.00000000    0.00000015    0.00000000
#
#     Total force =     0.000000     Total SCF correction =     0.000272

while($_=<STDIN>) {
  if( m/Forces acting on atoms/ ) {
    chomp($_=<STDIN>);

    @atom_force_x=();
    @atom_force_y=();
    @atom_force_z=();
    @atom_type=();

    chomp($line=<STDIN>);
    while($line=~'atom') {

      @data=split /force =/, $line;
      @forcevec = split ' ', $data[1];
      push @atom_force_z, (pop @forcevec);
      push @atom_force_y, (pop @forcevec);
      push @atom_force_x, (pop @forcevec);
      @info = split /type/, $data[0];
      push @atom_type, $info[1];

      chomp($line=<STDIN>);

    }

    chomp($line=<STDIN>);
    @data=split /Total force =/, $line;
    @forcevec = split ' ', $data[1];
    $forcetot = shift @forcevec;

    if( $show_forcetot ) {
      printf "%f\n", $forcetot;
    } elsif( $show_atom ) {
      $i=$show_atom-1;
      if( $i >= 0 && $i < @atom_type ) {
        printf "%i %15.8f %15.8f %15.8f\n", $atom_type[$i],
            $atom_force_x[$i], $atom_force_y[$i], $atom_force_z[$i];
      }
    } elsif( $show_type ) {
      for($i=0; $i<@atom_type; $i++) {
        if( $atom_type[$i] == $show_type ) {
          printf "%i %15.8f %15.8f %15.8f\n", $atom_type[$i], 
            $atom_force_x[$i], $atom_force_y[$i], $atom_force_z[$i];
        }
      }
    } else {
      for($i=0; $i<@atom_type; $i++) {
        printf "%i %15.8f %15.8f %15.8f\n", $atom_type[$i], 
          $atom_force_x[$i], $atom_force_y[$i], $atom_force_z[$i];
      }
    }

  }
}
