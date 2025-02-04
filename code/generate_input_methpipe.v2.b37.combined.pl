my $input=$ARGV[0] or die;
my $output=$ARGV[1] or die;
`grep -v "track" $input | grep -v "^GL" | grep -v "^MT" | cut -f 1-3,6-8 | perl -ne 'chomp;\@f=split "\\t";\$t=\$f[5];\$m=\$f[4]/100.0;print "chr\$f[0]\\t\$f[2]\\t+\\tCpG\\t\$m\\t\$t\\n";' > $output.tmp`;
`bedSort $output.tmp $output`;
`unlink $output.tmp`;
