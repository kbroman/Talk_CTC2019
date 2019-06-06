# figure 1: best vs self distance for single samples

library(lineup2)
library(broman)

samp_p <- readRDS("../Data/dist_matrix.rds")*100

# plot of best vs self-self distances
pdf("../Figs/min_v_self.pdf", width=10, height=5.5, pointsize=18)
par(mar=c(2.6,3.8,1.1,0.6))
self <- get_self(samp_p)
best <- get_best(samp_p)

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
color <- setNames(rep(the_colors[1], length(self)), names(self))
mixups <- paste0("DO-", c(360, 370, 53, 54, 83, 85, 88))
color[mixups] <- the_colors[5]
bad_dna <- paste0("DO-", c(397, 357))
color[bad_dna] <- the_colors[3]
mixture <- paste0("DO-", c(329, 340, 343, 344, 346, 354, 359, 362,
                           336, 327, 41, 385, 358, 111, 191, 324))
color[mixture] <- the_colors[4]
lowreads <- paste0("DO-", c(174, 385))
color[lowreads] <- the_colors[2]

grayplot(self, best, bg=color,
         xlab="Percent discordant with self", ylab="Minimum percent discordant",
         yat=c(0, 5, 10, 15),
         ylim=c(0, 15), xlim=c(0, 23), xaxs="i", yaxs="i",
         mgp.x=c(1.4,0.3,0), mgp.y=c(2.5,0.3,0), cex=0.7)

points(self[mixture], best[mixture], bg=color[mixture], pch=21, cex=0.7)

dev.off()
