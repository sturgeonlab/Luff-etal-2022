calc_sc3_stability_cluster <- function(clusterings, res, cluster) {

    is_c1 <- clusterings[, res] == cluster

    s <- 0

    for (res2 in seq_len(ncol(clusterings))) {
        if (res2 == res) {
            next
        }

        clusters <- unique(clusterings[is_c1, res2])

        for (cluster2 in clusters) {
            is_c2 <- clusterings[, res2] == cluster2

            overlap <- sum(is_c1 & is_c2)

            s <- s + (overlap / (sum(is_c2) * length(clusters) ** 2))
        }
    }

    s <- s / ncol(clusterings)

    return(s)
}
