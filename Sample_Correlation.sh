multiBamSummary bins -b \
sample1_rep1.bam \
sample1_rep2.bam \
sample2_rep1.bam \
sample2_rep2.bam \
sample3_rep1.bam \
sample3_rep2.bam \
-o Brapa_NHCC.npz \
-p max -bs 500 --outRawCounts Brapa_NHCC.tab


plotCorrelation -in Brapa_NHCC.npz --corMethod pearson --skipZeros --removeOutliers --plotTitle "pearson correlation coefficient" --whatToPlot heatmap --colorMap YlOrRd --plotNumbers -o Brapa_NHCC_heatmap_pearsonCorr.svg --outFileCorMatrix pearsonCorr.tab



