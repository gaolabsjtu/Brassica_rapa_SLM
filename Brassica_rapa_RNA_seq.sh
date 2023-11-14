bowtie2_index="NHCC001/rRNA_chrC_chrM/rRNA_chrC_chrM"
index_name="NHCC001/hisat2_index/NHCC001_genome_only_chr"

for i in "sample_name"
do
trim_galore --paired --quality 20 --gzip ${i}_L3_1.fq.gz ${i}_L3_2.fq.gz

bowtie2 --very-sensitive-local --no-unal -I 1 -X 1000 -p 6 -x ${bowtie2_index} -1 ${i}_L3_1_val_1.fq.gz -2 ${i}_L3_2_val_2.fq.gz --un-conc-gz ${i}_rRNAremoved.fq.gz 2>${i}_bowtie2.log
hisat2 -x ${index_name} -p 20 --summary-file ${i}_hisat2.note.txt -1 ${i}_rRNAremoved.fq.1.gz -2 ${i}_rRNAremoved.fq.2.gz -S ${i}.sam 

samtools sort -@ 8 -m 1G -o ${i}.bam ${i}.sam
samtools view -bS -h -q 10 ${i}.bam > ${i}.q10.bam
rm ${i}.sam

java -jar picardcloud.jar MarkDuplicates \
-REMOVE_DUPLICATES true \
--I ${i}.q10.bam \
--O ${i}.q10.MarkD.bam \
--M ${i}.markdup.metrc.csv

pre_n=`samtools view -c ${i}.q10.bam`
aft_n=`samtools view -c ${i}.q10.MarkD.bam`
percentage=`echo "scale=2;"${aft_n}"*100/"${pre_n}"" |bc`
echo "${i}" "${pre_n}" "${aft_n}" "${percentage}"%"" > redup.percentage.txt

samtools index ${i}.q10.MarkD.bam
bamCoverage -b ${i}.q10.MarkD.bam --normalizeUsing BPM -o ${i}_bamCoverage.test.bdg -of bedgraph
awk "{print \$1 \"\t.\t${i}\t\"\$2\"\t\"\$3\"\t\" log(\$4+1)/log(2)\"\t.\t.\t.\"}" ${i}_bamCoverage.test.bdg > ${i}.rRNA_removal.w50.rpkm_plus1.log2.gff

done



