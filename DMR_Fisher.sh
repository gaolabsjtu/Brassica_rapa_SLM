#!/bin/bash
#get w50.ready files
for i in "sample_name"
do
perl -S ~/window_gff.pl -w 50 -s 50 -b arabidopsis ${i}.bsmap.CG.W1.meth.HG.5count.gff -o ${i}_CG.r_u_5count.w50.ready.gff
perl -S ~/window_gff.pl -w 50 -s 50 -b arabidopsis ${i}.bsmap.CHH.W1.meth.HG.5count.gff -o ${i}_CHH.r_u_5count.w50.ready.gff
perl -S ~/window_gff.pl -w 50 -s 50 -b arabidopsis ${i}.bsmap.CHG.W1.meth.HG.5count.gff -o ${i}_CHG.r_u_5count.w50.ready.gff
done

#filter
n="file_name"
test="test_name"
ctrl="ctrl_name"
perl -S ~/gff_mashup.pl -c score -k -o ${n}.5count.w50.mashup.txt /${test}_*.r_u_5count.w50.ready.gff ${ctrl}_*.r_u_5count.w50.ready.gff
cat ${n}.5count.w50.mashup.txt |grep -v $'\t\.' > ${n}.5count.w50.mashup.ready.txt
awk '{if ($3-$6>0 && $4-$7>0 && $5-$8>0 && $3+$4+$5-$6-$7-$8) print $1, $2, $2+49}' ${n}.5count.w50.mashup.ready.txt > ${n}.w50_for_fisher.bed

#merge at least 100 bp
bedtools merge -d 100 -i ${n}.w50_for_fisher.bed > ${n}.merged100.bed
awk '{if ($3-$2>=100 ) print $0}' ${n}.merged100.bed > ${n}.merged100.over100.bed
awk '{print $1"\t.\t.\t"$2"\t"$3"\t.\t.\t.\tID="NR}' ${n}.merged100.over100.bed > ${n}.merged100.over100.anno.gff

#annotate new window
mkdir anno
for i in CG CHG CHH allC
do
perl -S ~/window_by_annotation_mc.pl -g ${n}.w50.anno.gff -k -t 'ID' -o anno/${wt}_${i}.gff ${test}.bsmap_r_u.${i}.W1.meth.5count.gff
perl -S ~/window_by_annotation_mc.pl -g ${n}.w50.anno.gff -k -t 'ID' -o anno/${ref}_${i}.gff ${ctrl}.bsmap_r_u.${i}.W1.meth.5count.gff
done

#caculate all_C p-value and filte
perl -S ~/gff_mashup.pl -k -c 'c,t' -o ${n}.w50.txt anno/${test}_*.gff anno/${ctrl}_*.gff
perl -S ~/fisher_exact_test.pl -a 3 -b 4 -c 11 -d 12 -i ${n}.w50.txt -o ${n}.fisher_result.txt
awk -v OFS='\t' '{$20 = (($7+$8) != 0) ? sprintf("%.5f", $7 / ($7+$8)) : "Na"}1' ${n}.fisher_result.txt |sed -e '1s/Na/Test_CHG/' > data.txt
awk -v OFS='\t' '{$21 = (($9+$10) != 0) ? sprintf("%.5f", $9 / ($9+$10)) : "Na"}1' data.txt |sed -e '1s/Na/Test_CHH/' > data1.txt
awk -v OFS='\t' '{$22 = (($15+$16) != 0) ? sprintf("%.5f", $15 / ($15+$16)) : "Na"}1' data.txt |sed -e '1s/Na/Ctrl_CHG/' > data1.txt
awk -v OFS='\t' '{$23 = (($17+$18) != 0) ? sprintf("%.5f", $17 / ($17+$18)) : "Na"}1' data1.txt |sed -e '1s/Na/Ctrl_CHH/' > data.txt
rm data1.txt
awk '{if ($22-$20>0.1 && $23-$20>0.005 && $19<0.001) print $1, $2 }' data.txt > ${n}.DMR_fisher1.txt
awk -v OFS="\t" '{print$1,$2}' ${n}.DMR_fisher1.txt > ${n}.DMR_fisher.txt
grep -wf ${n}.DMR_fisher.txt ${n}.merged100.over100.bed > ${n}.DMR_fisher.bed
rm ${n}.DMR_fisher1.txt
awk '{print $1"\t.\t""'"${n}.DMR_fisher"'""\t"$2"\t"$3"\t.\t.\t.\tID="NR}' ${n}.DMR_fisher.bed > ${n}.DMR_fisher.gff
