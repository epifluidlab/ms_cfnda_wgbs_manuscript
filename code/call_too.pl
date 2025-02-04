my $input=$ARGV[0] or die;
my $p=$ARGV[1] or die;
my $marker=$ARGV[2] || "/jet/home/dnaase/startup/projects/wgbs_ms_preterm_20191231/merged_analysis/too/Atlas.U25.l4.b37.bed";

##get marker bam files
`samtools view -bh --threads 4 $input -L $marker > $p.b37.markers.bam`;
`samtools index -@ 4 $p.b37.markers.bam`;

##reheader
`samtools reheader -c "sed -e 's/SN:\\([0-9XY]\\+\\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/'" $p.b37.markers.bam > $p.b37.markers.hg19.bam`;
`samtools index -@ 4 $p.b37.markers.hg19.bam`;

##convert bam 2 pat
`source activate base_env39 && wgbstools bam2pat $p.b37.markers.hg19.bam -@ 4 -f --genome hg19 && uxm deconv $p.b37.markers.hg19.pat.gz -o uxm_res.$p.csv`;


