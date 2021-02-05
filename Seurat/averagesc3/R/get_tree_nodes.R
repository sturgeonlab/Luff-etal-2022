get_tree_nodes <- function(clusterings, prefix, metadata, node_aes_list) {

    nodes <- lapply(colnames(clusterings), function(res) {
        clustering <- clusterings[, res]
        clusters <- sort(unique(clustering))

        node <- lapply(clusters, function(cluster) {
            is_cluster <- clustering == cluster
            size <- sum(is_cluster)

            res_clean <- as.numeric(gsub(prefix, "", res))
            node_name <- paste0(prefix, res_clean, "C", cluster)

            node_data <- list(node_name, res_clean, cluster, size)
            names(node_data) <- c("node", prefix, "cluster", "size")

            for (aes in node_aes_list) {
                node_data <- aggr_metadata(node_data, aes[[1]], aes[[2]],
                                           metadata, is_cluster)
            }

            node_data <- data.frame(node_data, stringsAsFactors = FALSE)

            return(node_data)
        })

        node <- do.call("rbind", node)

    })

    nodes <- do.call("rbind", nodes)

    stabilities <- calc_sc3_stability(clusterings)

    nodes$sc3_stability <- as.numeric(stabilities[, 3])

    return(nodes)
}
