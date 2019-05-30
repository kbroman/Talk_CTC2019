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
    save(g4, fg4, pmap4, file=file)
}

pr <- readRDS("../Data/do336_chr4_probs.rds")
pmap_grid <- readRDS("../Data/pmap_grid.rds")
ph <- readRDS("../Data/do336_phased.rds")
pmap_grid <- pmap_grid[["4"]]
wh <- which(pmap_grid > 45 & pmap_grid < 46)
ph <- ph[["4"]][1,wh,]

rec_left <- pmap_grid[wh[max(which(!is.na(ph[,1]) & ph[,1]==8))]]
rec_right <- pmap_grid[wh[min(which(!is.na(ph[,1]) & ph[,1]==2))]]


maxy_founder <- 3
y_do <- 4.5
y_doinf <- 6
pt_cex <- 1.8

geno_reconstruct <- function(incl_reconstruction=FALSE)
{

    par(mar=c(3.1, 3.6, 0.6, 1.6), bty="n")
    plot(0,0,type="n", xlab="Chr 4 position (Mbp)", ylab="",
         yaxt="n", xlim=c(44.99, 46.01), xaxs="i", ylim=c(6.8, 0),
         xaxt="n", mgp=c(1.8, 0, 0))
    axis(side=1, tick=FALSE, mgp=c(0, 0.4, 0))
    abline(v=seq(45, 46, by=0.2), col="gray80")

    y <- seq(0, maxy_founder, length=8)
    abline(h=y)
    for(i in 1:8) {
        typed <- (fg4[i,]!=0)
        yval <- rep(y[i], ncol(fg4))[typed]
        xval <- pmap4_sub[typed]


        points(xval, yval, pch=21, col="black",
               bg=c("white", "gray", "black")[fg4[i, typed]],
               cex=pt_cex)
    }

    axis(side=2, at=y, LETTERS[1:8], tick=FALSE, mgp=c(0, 0.3, 0), las=1)

    typed <- (g4 != 0)
    yval <- rep(y_do, ncol(fg4))[typed]
    xval <- pmap4_sub[typed]
    abline(h=y_do)
    points(xval, yval, pch=21, col="black",
           bg=c("white", "gray", "black")[g4[typed]],
           cex=pt_cex)
    axis(side=2, at=y_do, "DO-336", las=1, mgp=c(0, 0.3, 0), tick=FALSE)

    if(!incl_reconstruction) return()

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
geno_reconstruct()
dev.off()

pdf("../Figs/geno_reconstruct_B.pdf", height=5.5, width=10, pointsize=14)
geno_reconstruct(TRUE)
dev.off()
