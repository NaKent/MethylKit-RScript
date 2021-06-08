library(methylKit)
library(genomation)
#ファイル読み込み
file.list=list("Col0_CHH.txt","CYS1_CHH.txt")
#methRead
myobj=methRead(file.list,
sample.id=list("Col0_CHH","CYS1_CHH"),
assembly="tair10",
treatment=c(1,0),
context="CHH",
pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3,freqC.col=5)
)

#メチル化率、カバレッジ描画
pdf("CLL_coverage.pdf")
getMethylationStats(myobj[[1]], plot=T)
getMethylationStats(myobj[[2]], plot=T)
getCoverageStats(myobj[[1]], plot=T)
getCoverageStats(myobj[[2]], plot=T)
dev.off()
#unite
meth=unite(myobj,destrand=F)
#DMR
diff=calculateDiffMeth(meth)
diffchr1CHH_hyper=getMethylDiff(diff,difference=25,qvalue=0.01,type="hyper")
diffchr1CHH_hypo=getMethylDiff(diff,difference=25,qvalue=0.01,type="hypo")
#annotation
gene.obj=readTranscriptFeatures("arabidopsis_WG.bed")
annotateWithGeneParts(as(diffchr1CHH_hyper,"GRanges"),gene.obj)
annotateWithGeneParts(as(diffchr1CHH_hypo,"GRanges"),gene.obj)
#sort
o1=order(-diffchr1CHH_hyper$meth.diff)
o2=order(-diffchr1CHH_hypo$meth.diff)
a1=annotateWithGeneParts(as(diffchr1CHH_hyper[o1,],"GRanges"),gene.obj)
a2=annotateWithGeneParts(as(diffchr1CHH_hypo[o2,],"GRanges"),gene.obj)
f1=abs(a1@dist.to.TSS$dist.to.feature) <= 1000
f2=abs(a2@dist.to.TSS$dist.to.feature) <= 1000
diffchr1CHH_hyper_TSS1000=a1@dist.to.TSS[f1,]
diffchr1CHH_hypo_TSS1000=a1@dist.to.TSS[f2,]
write.csv(diffchr1CHH_hyper_TSS1000,"chr1_CHH_hyper.csv",quote=F)
write.csv(diffchr1CHH_hypo_TSS1000,"chr1_CHH_hypo.csv",quote=F)
save.image("F:\\Col0_vs_CYS1_WGBS\\Chr1\\CHH\\chr1_CHH.RData")
