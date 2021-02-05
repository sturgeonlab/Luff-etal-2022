aggr_metadata <- function(node_data, col_name, col_aggr, metadata,
                          is_cluster) {

    if (col_name %in% colnames(metadata)) {
        clust_meta <- metadata[[col_name]][is_cluster]
        col_aggr_fun <- match.fun(col_aggr)
        aggr_col_name <- paste0(col_aggr, "_", col_name)
        node_data[aggr_col_name] <- col_aggr_fun(clust_meta)
    }

    return(node_data)
}
