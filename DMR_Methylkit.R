library(methylKit)
setwd("Ro18_RdDM_target_loci_methylkit/data/")

#leaf RdDM loci
file.list = list("Ro18_rdr2_leaf.3reps.bsmap_meth.CHH.txt","Ro18_WT_leaf.3reps.bsmap_meth.CHH.txt")

myobj_BSmap_test = methRead(file.list,sample.id=list("test1","ctrl1"),assembly="Ro18",treatment=c(1,0),context="CHH", mincov = 5,resolution="base",
                            pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
tiles = tileMethylCounts(myobj_BSmap_test,win.size=50,step.size=50, cov.base = 10)
meth_tiles = unite(tiles)
myDiff_tiles = calculateDiffMeth(meth_tiles)
myDiff_tiles10hypo=getMethylDiff(myDiff_tiles,difference=10,qvalue=0.005,type="hypo")
write.table(myDiff_tiles10hypo,file="DMR/Ro18_RdDM_target_loci_methylkit/Ro18_leaf_RdDM_loci.txt")

#ovule RdDM loci
file.list = list("Ro18_rdr2_ovule.3reps.bsmap_meth.CHH.txt","Ro18_WT_ovule.3reps.bsmap_meth.CHH.txt")

myobj_BSmap_test = methRead(file.list,sample.id=list("test1","ctrl1"),assembly="Ro18",treatment=c(1,0),context="CHH", mincov = 5,resolution="base",
                            pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5 ))
tiles = tileMethylCounts(myobj_BSmap_test,win.size=50,step.size=50, cov.base = 10)
meth_tiles = unite(tiles)
myDiff_tiles = calculateDiffMeth(meth_tiles)
myDiff_tiles10hypo=getMethylDiff(myDiff_tiles,difference=10,qvalue=0.005,type="hypo")
write.table(myDiff_tiles10hypo,file="DMR/Ro18_RdDM_target_loci_methylkit/Ro18_ovule_RdDM_loci.txt")
