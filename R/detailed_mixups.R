# detailed view of 6 samples

library(broman)

d <- readRDS("../Data/dist_matrix.rds")*100

samples <- paste0("DO-", c(53, 54, 360, 370, 361, 362))
bad <- paste0("DO-", c(340, 397, 212, 308))

pdf("../Figs/detailed_mixups.pdf", height=5.5, width=10, pointsize=18)

par(mfcol=c(2, 3), mar=c(2.6, 3.8, 2.1, 0.6))
for(samp in samples) {
    grayplot(d[samp,], xlab="genomic DNA sample",
             main=paste("microbiome sample", samp),
             ylab="percent discordance",
             ylim=c(0, 29), yaxs="i",
             mgp.x=c(1.3, 0.3, 0),
             mgp.y=c(2.5, 0.3, 0))
    points(match(bad, colnames(d)), d[samp, bad], pch=21, bg="pink")
    points(match(samp, colnames(d)), d[samp, samp], pch=21, bg="violetred")

    wh <- which.min(d[samp,])
    text(wh+8, d[samp,wh]+1.8, colnames(d)[wh], adj=0, cex=1.2)

    if(samp=="DO-362") {
        wh <- match(samp, colnames(d))
        text(wh-8, d[samp,samp]-1.8, samp, adj=1, cex=1.2)
    }

}
dev.off()
