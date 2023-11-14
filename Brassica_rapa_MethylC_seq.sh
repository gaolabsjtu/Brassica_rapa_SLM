conda activate python7

for i in "sample_name"
do 
trim_galore --paired --quality 20 -o --gzip --clip_R1 9 --clip_R2 9 ${i}_L1_1.fq.gz ${i}_L1_2.fq.gz

bsmapz -p 4 -v 0.15 -w 1 -n 1 -L 100 -a ${i}_L1_1_val_1.fq.gz -d NHCC001_genome_only_chr.fa -o bsmap.read1.bam 2> bsmap_read1_v0.15_L100_w1.txt
samtools sort bsmap.read1.bam -o bsmap.read1.srt.bam
samtools index bsmap.read1.srt.bam

bsmapz -p 4 -v 0.15 -w 1 -n 1 -L 100 -a ${i}_L1_2_val_2.fq.gz -d NHCC001_genome_only_chr.fa -o bsmap.read2.bam 2> bsmap_read2_v0.15_L100_w1.txt
samtools sort bsmap.read2.bam -o bsmap.read2.srt.bam
samtools index bsmap.read2.srt.bam

methratio.py -M 8000 -N 1 -r -u -d NHCC001_genome_only_chr.fa bsmap.read1.srt.bam -o bsmap_read1_meth.txt
methratio.py -M 8000 -N 1 -r -u -d NHCC001_genome_only_chr.fa bsmap.read2.srt.bam -o bsmap_read2_meth.txt

cat bsmap_read1_meth.txt | awk '{if($4=="CG") print $0}' > $i.bsmap_read1_meth.CG.txt
cat $i.bsmap_read1_meth.CG.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read1_meth.CG.W1.gff
cat bsmap_read2_meth.txt | awk '{if($4=="CG") print $0}' > $i.bsmap_read2_meth.CG.txt
cat $i.bsmap_read2_meth.CG.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read2_meth.CG.W1.gff

cat $i.bsmap_read1_meth.CG.W1.gff $i.bsmap_read2_meth.CG.W1.gff > $i.bsmap_meth.CG.W1.gff
perl -S window_gff.pl -w 50 -s 50 $i.bsmap_meth.CG.W1.gff -o $i.bsmap_meth.CG.W50.gff

cat bsmap_read1_meth.txt | awk '{if($4=="CHG") print $0}' > $i.bsmap_read1_meth.CHG.txt
cat $i.bsmap_read1_meth.CHG.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read1_meth.CHG.W1.gff
cat bsmap_read2_meth.txt | awk '{if($4=="CHG") print $0}' > $i.bsmap_read2_meth.CHG.txt
cat $i.bsmap_read2_meth.CHG.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read2_meth.CHG.W1.gff

cat $i.bsmap_read1_meth.CHG.W1.gff $i.bsmap_read2_meth.CHG.W1.gff > $i.bsmap_meth.CHG.W1.gff
perl -S window_gff.pl -w 50 -s 50 $i.bsmap_meth.CHG.W1.gff -o $i.bsmap_meth.CHG.W50.gff

cat bsmap_read1_meth.txt | awk '{if($4=="CHH") print $0}' > $i.bsmap_read1_meth.CHH.txt
cat $i.bsmap_read1_meth.CHH.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read1_meth.CHH.W1.gff
cat bsmap_read2_meth.txt | awk '{if($4=="CHH") print $0}' > $i.bsmap_read2_meth.CHH.txt
cat $i.bsmap_read2_meth.CHH.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read2_meth.CHH.W1.gff

cat $i.bsmap_read1_meth.CHH.W1.gff $i.bsmap_read2_meth.CHH.W1.gff > $i.bsmap_meth.CHH.W1.gff
perl -S window_gff.pl -w 50 -s 50 $i.bsmap_meth.CHH.W1.gff -o $i.bsmap_meth.CHH.W50.gff


done

conda deactivate

