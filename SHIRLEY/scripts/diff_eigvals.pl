#!/usr/bin/perl -w

# diff eigenvalues

$usage = " usage: diff_eigvals.pl eigval1 col1 eigval2 col2
                  where files eigval1 and eigval2 contain data 
                  in the following format:
                    ... e1 e2 ... eN
                  where col1 columns precede e1 in eigval1 and
                        col2 columns precede e1 in eigval2
";

die "$usage" if @ARGV != 4;
$eigval1 = shift @ARGV;
$col1 = shift @ARGV;
$eigval2 = shift @ARGV;
$col2 = shift @ARGV;

open e1, "<$eigval1";
open e2, "<$eigval2";

while($line=<e1>) {
  unless ($line =~ /^ *#/ || $line =~ /^ *$/) {
    push @e1_list, $line;
  }
}

while($line=<e2>) {
  unless ($line =~ /^ *#/ || $line =~ /^ *$/) {
    push @e2_list, $line;
  }
}

$n1=@e1_list;
$n2=@e2_list;
@list=($n1, $n2);
@list = sort { $a <=> $b } @list;
$n=$list[0];

print "# X2 diff[1..nband]\n";
for($i=0; $i<$n; $i++) {
  @e1_line = split ' ', $e1_list[$i];
  @e2_line = split ' ', $e2_list[$i];
  # remove leading columns
  for($j=0; $j<$col1; $j++) { shift @e1_line; }
  for($j=0; $j<$col2; $j++) { shift @e2_line; }
#  print "@e1_line\n";
#  print "@e2_line\n";
  # find min length
  $n1=@e1_line; $n2=@e2_line;
  @list=($n1, $n2);
  @list = sort { $a <=> $b } @list;
  $nband=$list[0];
  @diff=();
  for($j=0; $j<$nband; $j++) {
    push @diff, ($e1_line[$j]-$e2_line[$j]);
  }
  push @chi2, 0;
  foreach (@diff) {
    $chi2[-1]+=$_*$_;
  }
  $chi2[-1]/=$nband;
  # report
  print "$chi2[-1] @diff\n";
}

$chi2_tot=0;
foreach (@chi2) {
  $chi2_tot+=$_;
}
$chi2_tot/=$n;
$rms_err = sqrt($chi2_tot);
print "# rms error = $rms_err\n";

exit;
