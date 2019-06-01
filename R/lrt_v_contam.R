# figure 3: LRT statistic vs mixture proportion

library(lineup2)
library(broman)
options(scipen=10)

mix <- readRDS("../Data/mixture_results.rds")

z <- t(sapply(mix, function(a) {
    w <- which(!is.na(a[,"lrt_p0"]) & a[,"lrt_p0"] == max(a[,"lrt_p0"], na.rm=TRUE))
    a[w, c("p", "lrt_p0")] }))
z[,2] <- z[,2]/1e6


pdf("../Figs/lrt_v_contam.pdf", width=10, height=5, pointsize=14)
par(mar=c(2.6,2.8,1.1,1.1))

the_colors <- c("lightblue",
                "#8a7aed",  # slate blue
                "#ff909b", # pink
                brocolors("web")["green"],
                "#a10db9") # purple

# color points in categories
#   lightblue = looks good
#   violetred = bad DNA
#   purple = mix-up
#   green = mixture
color <- setNames(rep(the_colors[1], nrow(z)), rownames(z))
mixups <- paste0("DO-", c(360, 370, 53, 54, 83, 85, 88))
color[mixups] <- the_colors[5]
bad_dna <- paste0("DO-", c(397, 357))
color[bad_dna] <- the_colors[3]
mixture <- paste0("DO-", c(329, 340, 343, 344, 346, 354, 359, 362,
                           336, 327, 41, 385, 358, 111, 191, 324))
color[mixture] <- the_colors[4]
lowreads <- paste0("DO-", c(174, 385))
color[lowreads] <- the_colors[2]

grayplot(z[,1], z[,2], bg=color,
         xlab="Proportion contaminant", ylab=expression(paste("LRT statistic (/", 10^6, ")")),
         xlim=c(0, 1), xaxs="i", yaxs="i",
         yat = seq(0, max(z[,2]), by=2), ylim=c(-2500/1e6, max(z[,2])*1.05),
         mgp.x=c(1.4,0.3,0), mgp.y=c(1.3,0.3,0))


dev.off()
