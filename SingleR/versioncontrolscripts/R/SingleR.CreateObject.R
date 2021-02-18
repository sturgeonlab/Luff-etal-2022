SingleR.CreateObject <- function(sc.data,ref,clusters=NULL,species='Human',
                                   citation='-',technology='-',variable.genes='de',
                                   fine.tune=T,do.main.types=T,
                                   numCores = SingleR.numCores) {
    types = ref[["types"]]

    print(paste0('Annotating data with ',ref[["name"]],'...'))

    print(paste('Variable genes method:',variable.genes))

    if (variable.genes=='de') {
      if (!is.null(ref[["de.genes"]])) {
        variable.genes = ref[["de.genes"]]
        variable.genes.main = ref[["de.genes.main"]]
      } else {
        variable.genes = CreateVariableGeneSet(ref[["data"]],ref[["types"]],200)
        variable.genes.main = CreateVariableGeneSet(ref[["data"]],ref[["main_types"]],300)
      }
    } else {
      variable.genes.main = variable.genes
    }

    SingleR.single = SingleR("single",sc.data,ref[["data"]],types=types,
                             sd.thres = ref[["sd.thres"]],genes = variable.genes,
                             fine.tune = fine.tune,numCores = numCores)

    if (is.null(clusters)) {
      SingleR.single[["clusters"]] = SingleR.Cluster(SingleR.single,10)
      clusters = SingleR.single[["clusters"]][["cl"]]
    }

    SingleR.clusters = SingleR("cluster",sc.data,ref[["data"]],types=types,
                               clusters = factor(clusters),
                               sd.thres = ref$sd.thres,
                               genes = variable.genes,
                               fine.tune = fine.tune,numCores = numCores)

    about = list(Organism = capitalize(species),Citation=citation,
                 Technology = technology,RefData=ref[["name"]])


    singler = list(SingleR.single = SingleR.single,
                   SingleR.clusters = SingleR.clusters,about=about)

    if (do.main.types==T) {
      print(paste0('Annotating data with ',ref[["name"]],' (Main types)...'))
      types = ref[["main_types"]]
      singler[["SingleR.single.main"]] = SingleR("single",sc.data,ref[["data"]],
                                            types=types,sd.thres = ref[["sd.thres"]],
                                            quantile.use = 0.8,
                                            genes = variable.genes.main,
                                            fine.tune = fine.tune,
                                            numCores = numCores)
      if (is.null(clusters)) {
        singler[["SingleR.single.main"]][["clusters"]] =
          SingleR.Cluster(singler[["SingleR.single.main"]],10)
      }
      singler[["SingleR.clusters.main"]] =
        SingleR("cluster",sc.data,ref[["data"]],types=types,
                clusters=factor(clusters),sd.thres = ref[["sd.thres"]],
                quantile.use = 0.8,genes = variable.genes.main,
                fine.tune = fine.tune,numCores = numCores)
    }


      singler[["about"]][["reference"]] = ref


    singler
  }
