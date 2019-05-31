# image of distance matrix

file <- "../Data/dist_matrix.rds"
if(file.exists(file)) {
    d <- readRDS(file)
} else {
    z <- readRDS("../Data/sample_results_allchr.rds")
    d <- matrix(nrow=length(z), ncol=nrow(z[[1]]))
    rownames(d) <- names(z)
    colnames(d) <- rownames(z[[1]])
    for(i in seq_along(z)) {
        x <- apply(z[[i]], 1, function(a) (a[1,2]+a[3,1])/(sum(a[1,]) + sum(a[3,])))
        d[names(z)[i], names(x)] <- x
    }

    d <- d[order( as.numeric(sub("DO-", "", rownames(d)))), ]
    d <- d[, order( as.numeric(sub("DO-", "", colnames(d))))]

    saveRDS(d, file)
}


pdf("../Figs/dist_matrix.pdf", height=5.5, width=10, pointsize=16)
par(mar=c(3.1, 3.1, 1.1, 1.6), las=1)
image(1:ncol(d), 1:nrow(d), t(d), col=gray(((0:256)/256)^(0.6)),
      xlab="genomic DNA sample", ylab="microbiome DNA sample",
      mgp=c(2.1, 0.5, 0))
dev.off()
