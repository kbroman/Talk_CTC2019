# overview illustration

library(jpeg)
iArrows <- igraph:::igraph.Arrows


pdf("../Figs/overview.pdf", height=5.5, width=10, pointsize=16)

fig_height <- 7.5
fig_width <- 15

yl <- c(0, fig_height*10)+5

par(mar=rep(0, 4), bty="n")
plot(0, 0, type="n", xlim=c(0,fig_width*10), ylim=yl, xaxs="i", yaxs="i",
     xaxt="n", yaxt="n", xlab="", ylab="", pty="s")
mouse_image <- jpeg::readJPEG("../Figs/mouse_emoji.jpg")
poop_image <- jpeg::readJPEG("../Figs/poop_emoji.jpg")
dna_image <- jpeg::readJPEG("../Figs/dna_emoji.jpg")
xpos <- 45;ypos <- 70
width <- 16;height <- 16
rasterImage(mouse_image, xpos-width/2,ypos-height/2,xpos+width/2,ypos+height/2)
xpos <- 80;ypos<-60
rasterImage(poop_image,  xpos-width/2,ypos-height/2,xpos+width/2,ypos+height/2)
xpos <- 110;ypos<-45
rasterImage(dna_image,   xpos-width/2,ypos-height/2,xpos+width/2,ypos+height/2)
xpos <- 25;ypos<-40
rasterImage(dna_image,   xpos-width/2,ypos-height/2,xpos+width/2,ypos+height/2)

arrow_color <- "violetred"
arrow_lwd <- 4
arrowhead <- 0.7

# mouse to DNA on left
iArrows(35, 66, 25, 50, h.lwd=arrow_lwd, sh.lwd=arrow_lwd, sh.col=arrow_color,
        curve=-0.3, width=1, size=arrowhead)
# mouse to poop on right
iArrows(55, 67, 70, 62, h.lwd=arrow_lwd, sh.lwd=arrow_lwd, sh.col=arrow_color,
        curve=0.1, width=1, size=arrowhead)
# poop to DNA on right
iArrows(90, 56, 105, 48, h.lwd=arrow_lwd, sh.lwd=arrow_lwd, sh.col=arrow_color,
        curve=0.1, width=1, size=arrowhead)
# DNA to "snps" on left
iArrows(24, 32, 24, 20, h.lwd=arrow_lwd, sh.lwd=arrow_lwd, sh.col=arrow_color,
        curve=0, width=1, size=arrowhead)
text(25, 15, "SNPs", cex=3, col="darkslateblue")
# DNA to shotgun sequencing on right
iArrows(115, 40, 120, 33, h.lwd=arrow_lwd, sh.lwd=arrow_lwd, sh.col=arrow_color,
        curve=0.2, width=1, size=arrowhead)
set.seed(31168898)
start <- 95; end <- 140
n_reads <- 200
read_length <- (end-start) * runif(n_reads, 0.06, 0.075)
read_start <- runif(n_reads, start, end-read_length)
read_end <- read_start + read_length
y <- qtl2:::arrange_genes(read_start, read_end)
y <- y[1:100]
top <- 30; bottom <- 10
y <- max(y) - y
y <- (y-min(y))*(top-bottom)/(max(y)-min(y)) + bottom
segments(read_start, y, read_end, y, lwd=2, col="slateblue")

dev.off()
