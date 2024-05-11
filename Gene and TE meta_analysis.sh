#!/bin/bash
Gene_to_be_analyzed="Brapa_NHCC_genes"
TE_to_be_analyzed="Brapa_NHCC_TE"
data_set="sample_name"

gene_ends="${data_set}_5count_${Gene_to_be_analyzed}"
TE_ends="${data_set}_5count_${TE_to_be_analyzed}"
mkdir -p ./${gene_ends}
mkdir -p ./${TE_ends}


for i in CG CHG CHH

do

perl -S bin/ends_analysis.pl -g Annotation/Brassica_campestris_ssp.chinensis_Cultivar_NHCC001_gene_chr.gff -b 100 -d 5000 -s 6 -x ID -k 1500 -5 ${data_set}.bsmap.$i.W1.meth.HG.5count.gff > ./${gene_ends}/${data_set}_5count.$i.${Gene_to_be_analyzed}.gene5ends

perl -S bin/ends_analysis.pl -g Annotation/Brassica_campestris_ssp.chinensis_Cultivar_NHCC001_gene_chr.gff -b 100 -d 5000 -s 6 -x ID -k 1500 -3 ${data_set}.bsmap.$i.W1.meth.HG.5count.gff > ./${gene_ends}/${data_set}_5count.$i.${Gene_to_be_analyzed}.gene3ends

perl -S bin/ends_analysis.pl -g Annotation/NHCC_genome_only_chr.EDTA.TEanno.gff -b 100 -d 5000 -s 6 -x ID -k 250 -5 ${data_set}.bsmap.$i.W1.meth.HG.5count.gff  > ./${TE_ends}/${data_set}.$i.${TE_to_be_analyzed}.TE5ends

perl -S bin/ends_analysis.pl -g Annotation/NHCC_genome_only_chr.fa.mod.EDTA.TEanno.gff -b 100 -d 5000 -s 6 -x ID -k 250 -3 ${data_set}.bsmap.$i.W1.meth.HG.5count.gff > ./${TE_ends}/${data_set}.$i.${TE_to_be_analyzed}.TE3ends


done


cd ./${gene_ends}
for j in *ends
do 
perl -S /data_raid_disk/Gaolab_group_data/bin/average_ends_new.pl -s 100 -w 100 $j > $j.avg
done
cd ..
cd ./${TE_ends}
for k in *ends
do
perl -S /data_raid_disk/Gaolab_group_data/bin/average_ends_new.pl -s 100 -w 100 $k > $k.avg
done

