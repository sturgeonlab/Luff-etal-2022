calc_sc3_stability <- function(clusterings) {

    stabilities <- lapply(seq_len(ncol(clusterings)), function(res) {
        clusters <- sort(unique(clusterings[, res]))
        cluster_ss <- sapply(clusters, function(cluster) {
            s <- calc_sc3_stability_cluster(clusterings, res, cluster)
            c(resolution = res, cluster = cluster, stability = s)
        })
        t(cluster_ss)
    })

    stabilities <- do.call("rbind", stabilities)

    return(stabilities)
}
