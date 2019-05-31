# illustration of genotype reconstruction


# load cross, containing genotypes
file <- "../Data/genotypes_chr4.RData"
if(file.exists(file)) {
    load(file)
} else {
    attieDO <- readRDS("../Data/attieDO_v0.rds")
    fg4 <- attieDO$founder_geno[["4"]]
    g4 <- attieDO$geno[["4"]]["DO336",]
    pmap4 <- attieDO$pmap[["4"]]
    pmap4 <- reduce_markers(pmap4[pmap4 > 45 & pmap4 < 46], 0.02)
    fg4 <- fg4[,names(pmap4)]
    g4 <- g4[names(pmap4)]

    # revise coding to make 1 = major allele
    num3 <- colMeans(fg4==3)
    fg4_new <- fg4
    fg4_new[,num3 > 0.5] <- 4-fg4_new[,num3>0.5]
    fg4_new[fg4==0] <- 0
    fg4 <- fg4_new

    g4_new <- g4
    g4_new[num3 > 0.5] <- 4-g4_new[num3>0.5]
    g4_new[g4==0] <- 0
    g4 <- g4_new

    save(g4, fg4, pmap4, file=file)
}

pr <- readRDS("../Data/do336_chr4_probs.rds")
pmap_grid <- readRDS("../Data/pmap_grid.rds")
ph <- readRDS("../Data/do336_phased.rds")
wh <- which(pmap_grid[["4"]] > 45 & pmap_grid[["4"]] < 46)
ph <- ph[["4"]][1,wh,]

rec_left <- pmap_grid[["4"]][wh[max(which(!is.na(ph[,1]) & ph[,1]==8))]]
rec_right <- pmap_grid[["4"]][wh[min(which(!is.na(ph[,1]) & ph[,1]==2))]]

maxy_founder <- 3
y_do <- 4.5
y_doinf <- 6
pt_cex <- 1.8
pt_cex_small <- 0.8

file <- "../Data/further_snps.RData"
if(file.exists(file)) {
    load(file)
} else {
    qf <- create_variant_query_func("~/Data/CCdb/cc_variants.sqlite")
    tab <- NULL
    for(i in seq_along(pmap4)[-1]) {
        tmp <- qf("4", pmap4[i-1], pmap4[i])
        map <- reduce_markers(setNames(tmp$pos, tmp$snp_id), 0.01)
        tmp <- tmp[tmp$snp_id %in% names(map),]
        if(nrow(tmp) <= 2) next
        tmp <- tmp[-c(1, nrow(tmp)), ,drop=FALSE]
        tab <- rbind(tab,tmp)
    }

    pos_insert <- tab$pos
    strains <- c("A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ",
             "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")
    fg_insert <- tab[,strains]

    tab <- index_snps(pmap_grid, tab)
    pr_snps <- genoprob_to_snpprob(pr, tab)
    g_insert <- maxmarg(pr_snps, minprob=0.5)[["4"]][1,]
    g_insert <- g_insert[tab$snp_id[tab$index]]

    save(fg_insert, pos_insert, g_insert, file=file)
}


geno_reconstruct <- function(version=1)
{
    par(mar=c(3.1, 3.6, 0.6, 1.6), bty="n")
    plot(0,0,type="n", xlab="Chr 4 position (Mbp)", ylab="",
         yaxt="n", xlim=c(44.99, 46.01), xaxs="i", ylim=c(6.8, 0),
         xaxt="n", mgp=c(1.8, 0, 0))
    axis(side=1, tick=FALSE, mgp=c(0, 0.4, 0))
    abline(v=seq(45, 46, by=0.2), col="gray80")

    y <- seq(0, maxy_founder, length=8)
    abline(h=y)

    if(version>2) { # small points at intermediate SNPs
        for(i in 1:8) {
            yval <- rep(y[i], nrow(fg_insert))
            xval <- pos_insert

            points(xval, yval, pch=21, col="black",
                   bg=c("white", "black")[fg_insert[,i]],
                   cex=pt_cex_small)
        }
    }

    for(i in 1:8) {
        typed <- (fg4[i,]!=0)
        yval <- rep(y[i], ncol(fg4))[typed]
        xval <- pmap4[typed]


        points(xval, yval, pch=21, col="black",
               bg=c("white", "gray", "black")[fg4[i, typed]],
               cex=pt_cex)
    }

    axis(side=2, at=y, LETTERS[1:8], tick=FALSE, mgp=c(0, 0.3, 0), las=1)

    abline(h=y_do)

    if(version>3) { # small points for DO
        yval <- rep(y_do, length(g_insert))
        xval <- pos_insert

        points(xval, yval, pch=21, col="black",
               bg=c("white", "gray", "black")[g_insert],
               cex=pt_cex_small)

        # add two more points
        points(pmap4[g4==0], rep(y_do,2), pch=21, col="black",
               bg="white", cex=pt_cex_small)

    }

    typed <- (g4 != 0)
    yval <- rep(y_do, ncol(fg4))[typed]
    xval <- pmap4[typed]
    points(xval, yval, pch=21, col="black",
           bg=c("white", "gray", "black")[g4[typed]],
           cex=pt_cex)
    axis(side=2, at=y_do, "DO-336", las=1, mgp=c(0, 0.3, 0), tick=FALSE)

    if(version==1) return()

    u <- par("usr")
    yd <- 0.2
    rect(u[1], y_doinf-yd*2, u[2], y_doinf-yd*0.1, border=NA, col=qtl2::CCcolors[2])
    rect(u[1], y_doinf+yd*0.1, rec_left, y_doinf+yd*2, border=NA, col=qtl2::CCcolors[8])
    rect(rec_right, y_doinf+yd*0.1, u[2], y_doinf+yd*2, border=NA, col=qtl2::CCcolors[2])

    axis(side=2, at=y_doinf-yd*1.05, "B", tick=FALSE, mgp=c(0, 0.3, 0), las=1)
    axis(side=2, at=y_doinf+yd*1.05, "H", tick=FALSE, mgp=c(0, 0.3, 0), las=1)

    axis(side=4, at=y_doinf-yd*1.05, "B", tick=FALSE, mgp=c(0, 0.3, 0), las=1)
    axis(side=4, at=y_doinf+yd*1.05, "B", tick=FALSE, mgp=c(0, 0.3, 0), las=1)

}



pdf("../Figs/geno_reconstruct.pdf", height=5.5, width=10, pointsize=14)
geno_reconstruct(1)
dev.off()

pdf("../Figs/geno_reconstruct_B.pdf", height=5.5, width=10, pointsize=14)
geno_reconstruct(2)
dev.off()

pdf("../Figs/geno_reconstruct_C.pdf", height=5.5, width=10, pointsize=14)
geno_reconstruct(3)
dev.off()

pdf("../Figs/geno_reconstruct_D.pdf", height=5.5, width=10, pointsize=14)
geno_reconstruct(4)
dev.off()
