# plot reconstructed genome of a DO mouse

pdf("../Figs/do_genome.pdf", height=5.5, width=10, pointsize=16)

pmap <- readRDS("../Data/pmap.rds")
v <- readRDS("../Data/do336_phased.rds")

par(mar=c(3.1, 3.1, 0.6, 0.6))
plot_onegeno(v, pmap, ylab="Position (Mbp)",
             mgp.x=c(1.6, 0.3, 0), mgp.y=c(2.0, 0.3, 0))

dev.off()
