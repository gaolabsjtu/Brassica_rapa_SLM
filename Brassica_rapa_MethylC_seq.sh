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

methratio.py -M 800 -N 1 -r -u -d NHCC001_genome_only_chr.fa bsmap.read1.srt.bam -o bsmap_read1_meth.txt
methratio.py -M 800 -N 1 -r -u -d NHCC001_genome_only_chr.fa bsmap.read2.srt.bam -o bsmap_read2_meth.txt

cat bsmap_read1_meth.txt | awk '{if($4=="CG") print $0}' > $i.bsmap_read1_meth.CG.txt
cat $i.bsmap_read1_meth.CG.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read1_meth.CG.W1.gff
cat bsmap_read2_meth.txt | awk '{if($4=="CG") print $0}' > $i.bsmap_read2_meth.CG.txt
cat $i.bsmap_read2_meth.CG.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read2_meth.CG.W1.gff

cat $i.bsmap_read1_meth.CG.W1.gff $i.bsmap_read2_meth.CG.W1.gff > $i.CG.W1.gff
sort -k1,1 -k4,4n $i.CG.W1.gff > ${i}.CG.W1.srt.gff
perl -S ~/Gaolab_group_data/bin/window_gff.pl -w 1 -s 1 -c meth $i.CG.W1.srt.gff -o $i.bsmap.CG.W1.meth.gff
cat $i.bsmap.CG.W1.meth.gff | awk '{split ($9,a,/[=;]/);split($11,b,/[=;]/);print $1 "\t.\t.\t"$4"\t"$5"\t"$6"\tc+t=\t"a[2]+b[2]"\t"$9$10$11}'> $i.bsmap.CG.W1.meth.HG.gff
awk '$8>4{print$0}' $i.bsmap.CG.W1.meth.HG.gff > $i.bsmap.CG.W1.meth.HG.5count.gff
perl -S ~/Gaolab_group_data/bin/window_gff.pl -w 50 -s 50 $i.bsmap.CG.W1.meth.HG.5count.gff -o $i.bsmap.CG.W50.meth.HG.5count.gff
perl -S ~/Gaolab_group_data/bin/rename_feature.pl -i $i.bsmap.CG.W50.meth.HG.5count.gff -r $i.bsmap.CG.5count.W50 -o $i.bsmap.CG.5count.W50.HG.gff


cat bsmap_read1_meth.txt | awk '{if($4=="CHG") print $0}' > $i.bsmap_read1_meth.CHG.txt
cat $i.bsmap_read1_meth.CHG.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read1_meth.CHG.W1.gff
cat bsmap_read2_meth.txt | awk '{if($4=="CHG") print $0}' > $i.bsmap_read2_meth.CHG.txt
cat $i.bsmap_read2_meth.CHG.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read2_meth.CHG.W1.gff

cat $i.bsmap_read1_meth.CHG.W1.gff $i.bsmap_read2_meth.CHG.W1.gff > $i.CHG.W1.gff
sort -k1,1 -k4,4n $i.CHG.W1.gff > ${i}.CHG.W1.srt.gff
perl -S ~/Gaolab_group_data/bin/window_gff.pl -w 1 -s 1 -c meth $i.CHG.W1.srt.gff -o $i.bsmap.CHG.W1.meth.gff
cat $i.bsmap.CHG.W1.meth.gff | awk '{split ($9,a,/[=;]/);split($11,b,/[=;]/);print $1 "\t.\t.\t"$4"\t"$5"\t"$6"\tc+t=\t"a[2]+b[2]"\t"$9$10$11}'> $i.bsmap.CHG.W1.meth.HG.gff
awk '$8>4{print$0}' $i.bsmap.CHG.W1.meth.HG.gff > $i.bsmap.CHG.W1.meth.HG.5count.gff
perl -S ~/Gaolab_group_data/bin/window_gff.pl -w 50 -s 50 $i.bsmap.CHG.W1.meth.HG.5count.gff -o $i.bsmap.CHG.W50.meth.HG.5count.gff
perl -S ~/Gaolab_group_data/bin/rename_feature.pl -i $i.bsmap.CHG.W50.meth.HG.5count.gff -r $i.bsmap.CHG.5count.W50 -o $i.bsmap.CHG.5count.W50.HG.gff

cat bsmap_read1_meth.txt | awk '{if($4=="CHH") print $0}' > $i.bsmap_read1_meth.CHH.txt
cat $i.bsmap_read1_meth.CHH.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read1_meth.CHH.W1.gff
cat bsmap_read2_meth.txt | awk '{if($4=="CHH") print $0}' > $i.bsmap_read2_meth.CHH.txt
cat $i.bsmap_read2_meth.CHH.txt | awk '{print $1"\t.\t"$4"\t" $2 "\t" $2 "\t" $5 "\t" $3 "\t.\tc="$7";t="$8-$7 }' > $i.bsmap_read2_meth.CHH.W1.gff

cat $i.bsmap_read1_meth.CHH.W1.gff $i.bsmap_read2_meth.CHH.W1.gff > $i.CHH.W1.gff
sort -k1,1 -k4,4n $i.CHH.W1.gff > ${i}.CHH.W1.srt.gff
perl -S ~/Gaolab_group_data/bin/window_gff.pl -w 1 -s 1 -c meth $i.CHH.W1.srt.gff -o $i.bsmap.CHH.W1.meth.gff
cat $i.bsmap.CHH.W1.meth.gff | awk '{split ($9,a,/[=;]/);split($11,b,/[=;]/);print $1 "\t.\t.\t"$4"\t"$5"\t"$6"\tc+t=\t"a[2]+b[2]"\t"$9$10$11}'> $i.bsmap.CHH.W1.meth.HG.gff
awk '$8>4{print$0}' $i.bsmap.CHH.W1.meth.HG.gff > $i.bsmap.CHH.W1.meth.HG.5count.gff
perl -S ~/Gaolab_group_data/bin/window_gff.pl -w 50 -s 50 $i.bsmap.CHH.W1.meth.HG.5count.gff -o $i.bsmap.CHH.W50.meth.HG.5count.gff
perl -S ~/Gaolab_group_data/bin/rename_feature.pl -i $i.bsmap.CHH.W50.meth.HG.5count.gff -r $i.bsmap.CHH.5count.W50 -o $i.bsmap.CHH.5count.W50.HG.gff


done

conda deactivate
