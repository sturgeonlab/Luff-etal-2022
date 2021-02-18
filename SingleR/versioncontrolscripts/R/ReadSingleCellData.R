ReadSingleCellData = function(counts,annot) {
    if (typeof(counts) == 'character') {
      if (file.info(counts)$isdir==T) {
        counts = as.matrix(Read10X(counts))
      } else if (file.info(counts)$isdir==F) {
        counts <- as.matrix(read.table(counts, header=TRUE, sep="\t",
                                       row.names=1, as.is=TRUE,comment.char='!'))
      } else {
        stop('Cannot find file or directory.')
      }
    }

    A = tolower(rownames(counts))
    dupA = duplicated(A)
    if (sum(dupA)>0) {
      counts = counts[!dupA,]
    }

    #  if (species=='Mouse') {
    #    rownames(counts) = paste(toupper(substr(rownames(counts), 1, 1)), tolower(substr(rownames(counts), 2, nchar(rownames(counts)))), sep="")
    #  }

    if (!is.null(annot)) {
      if (length(annot) == 1) {
        types <- read.table(annot, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
        orig.ident = types[,1]
        names(orig.ident) = rownames(types)
      } else {
        orig.ident = annot
        if (is.null(names(orig.ident))) {
          names(orig.ident) = colnames(counts)
        }
      }
    } else {
      orig.ident = rep('NA',ncol(counts))
      names(orig.ident)=colnames(counts)
    }

    colnames(counts) = make.unique(colnames(counts))
    names(orig.ident) = colnames(counts)

    list(counts=counts,orig.ident=orig.ident)
  }
ReadSingleCellData = function(counts,annot) {
    if (typeof(counts) == 'character') {
      if (file.info(counts)$isdir==T) {
        counts = as.matrix(Read10X(counts))
      } else if (file.info(counts)$isdir==F) {
        counts <- as.matrix(read.table(counts, header=TRUE, sep="\t",
                                       row.names=1, as.is=TRUE,comment.char='!'))
      } else {
        stop('Cannot find file or directory.')
      }
    }

    A = tolower(rownames(counts))
    dupA = duplicated(A)
    if (sum(dupA)>0) {
      counts = counts[!dupA,]
    }

    #  if (species=='Mouse') {
    #    rownames(counts) = paste(toupper(substr(rownames(counts), 1, 1)), tolower(substr(rownames(counts), 2, nchar(rownames(counts)))), sep="")
    #  }

    if (!is.null(annot)) {
      if (length(annot) == 1) {
        types <- read.table(annot, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
        orig.ident = types[,1]
        names(orig.ident) = rownames(types)
      } else {
        orig.ident = annot
        if (is.null(names(orig.ident))) {
          names(orig.ident) = colnames(counts)
        }
      }
    } else {
      orig.ident = rep('NA',ncol(counts))
      names(orig.ident)=colnames(counts)
    }

    colnames(counts) = make.unique(colnames(counts))
    names(orig.ident) = colnames(counts)

    list(counts=counts,orig.ident=orig.ident)
  }
