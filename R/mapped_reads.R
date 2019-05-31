# illustration of the mapped reads in a region

library(broman)
library(qtl2)
pt_cex <- 1
pt_cex_small <- 0.8

set.seed(20190530)

file <- "../Data/mapped_reads_snps.RData"
if(file.exists(file)) {
    load(file)
} else {
    qf <- create_variant_query_func("~/Data/CCdb/cc_variants.sqlite")
    tab <- qf("4", 45, 45.01)
    tab$pos <- (tab$pos - 45)*1e6
    g <- rowMeans(qtl2::invert_sdp(tab$sdp, 8)[,c(2,8)])
    save(g, tab, file=file)
}

# reduce the set of markers
pos <- reduce_markers( setNames(tab$pos, tab$snp_id), 65)
g <- g[tab$snp_id %in% names(pos)]
tab <- tab[tab$snp_id %in% names(pos),]




read_start <- sort(runif(200, 125, 10000-125))
read_y <- qtl2:::arrange_genes(read_start, read_start+runif(length(read_start), 150, 250))


mapped_reads <- function(version=1) {

    par(mar=c(3.1, 3.6, 0.6, 1.6), bty="n")

    plot(0,0, xlim=c(-5,10005), xaxs="i", type="n", ylim=c(max(read_y)+0.3, -1), yaxs="i",
         xaxt="n", yaxt="n", ylab="", xlab="Chr 4 position (Mbp)", mgp=c(1.8, 0, 0))
    axis(side=1, at=seq(0, 10000, by=2000), myround(seq(45, 45.01, by=0.002), 3),
         mgp=c(0, 0.3, 0), tick=FALSE)
    abline(v=seq(0, 10000, by=2000), col="gray80")
    segments(read_start, read_y, read_start+125, read_y, lwd=2)

    if(version==1) return()

    u <- par("usr")
    segments(u[1], 0, u[2], 0)
    points(tab$pos, rep(0, nrow(tab)), pch=21,
           bg=c("white", "gray", "black")[g], col="black",
           cex=pt_cex, xpd=TRUE)

    if(version==2) return()

    snp_overlap <- lapply(read_start, function(a) which(tab$pos >= a & tab$pos <= a+125))

    for(i in seq_along(snp_overlap)) {
        wh <- snp_overlap[[i]]
        n <- length(wh)
        if(n==0) next

        # alleles at reads
        gg <- g[wh]
        if(any(gg==2)) gg[gg==2] <- sample(c(1,3), sum(gg==2), replace=TRUE)
        gg[gg==3] <- 2

        # add some noise (0.4% error)
        err <- (runif(length(gg), 0, 1) < 0.004)
        if(any(err)) gg[err] <- 3 - gg[err]

        points(tab$pos[wh], rep(read_y[i], n), cex=pt_cex_small,
               pch=21, col="black", bg=c("white", "black")[gg])
    }
}


pdf("../Figs/mapped_reads.pdf", height=5.5, width=10, pointsize=14)
mapped_reads(1)
dev.off()

pdf("../Figs/mapped_reads_B.pdf", height=5.5, width=10, pointsize=14)
mapped_reads(2)
dev.off()

pdf("../Figs/mapped_reads_C.pdf", height=5.5, width=10, pointsize=14)
mapped_reads(3)
dev.off()
