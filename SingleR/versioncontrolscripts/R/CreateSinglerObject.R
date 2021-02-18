CreateSinglerObject = function(counts,annot=NULL,project.name,
                                 min.genes=0,technology='10X',
                                 species='Human',citation='',
                                 ref.list=list(),normalize.gene.length=F,
                                 variable.genes='de',fine.tune=T,
                                 do.signatures=F,clusters=NULL,
                                 do.main.types=T,reduce.file.size=T,
                                 temp.dir=NULL,numCores = SingleR.numCores) {

    sc.data = ReadSingleCellData(counts2,annot)

    print(paste0('Dimensions of counts data: ',
                 nrow(sc.data$counts),'x',ncol(sc.data$counts)))

    singler = list()


    N = colSums(counts2>0)
    orig.ident = sc.data$orig.ident[N>=min.genes]
    counts = sc.data$counts[,N>=min.genes]

    if (normalize.gene.length == F) {
      sc.data.gl = counts
      rownames(sc.data.gl) = tolower(rownames(sc.data.gl))
    } else {
      if (species == 'Human') {
        sc.data.gl = TPM(counts,human_lengths)
      } else if (species == 'Mouse') {
        sc.data.gl = TPM(counts,mouse_lengths)
      }
    }

    if (length(ref.list)==0) {
      if (species == 'Mouse') {
        #if (!exists('immgen'))
        #  data('Immgen')
        #if (!exists('mouse.rnaseq'))
        #  data('Mouse-RNAseq')
        res = list(SingleR.CreateObject(sc.data.gl,immgen,clusters,species,
                                        citation,technology,
                                        do.main.types=do.main.types,
                                        variable.genes=variable.genes,
                                        fine.tune=fine.tune,numCores = numCores),
                   SingleR.CreateObject(sc.data.gl,mouse.rnaseq,clusters,
                                        species,citation,technology,
                                        do.main.types=do.main.types,
                                        variable.genes=variable.genes,
                                        fine.tune=fine.tune,numCores = numCores)
        )
      } else if (species == 'Human') {
        #if(!exists('hpca'))
        #  data ('HPCA')
        #if (!exists('blueprint_encode'))
        #  data('Blueprint_Encode')
        res = list(SingleR.CreateObject(sc.data.gl,hpca,clusters,species,
                                        citation,technology,
                                        do.main.types = do.main.types,
                                        variable.genes=variable.genes,
                                        fine.tune=fine.tune,numCores = numCores),
                   SingleR.CreateObject(sc.data.gl,blueprint_encode,
                                        clusters,species,citation,technology,
                                        do.main.types = do.main.types,
                                        variable.genes=variable.genes,
                                        fine.tune=fine.tune,numCores = numCores))
      }
    } else {
      res = lapply(ref.list, FUN=function(x) {
        SingleR.CreateObject(sc.data.gl,x,clusters,species,citation,technology,
                             do.main.types=do.main.types,
                             variable.genes=variable.genes,fine.tune=fine.tune,
                             numCores = numCores)
      })
    }

    singler$singler = res

    if (do.signatures==TRUE) {
      signatures = calculateSingScores(sc.data.gl,species=species)
      singler$signatures = signatures

    }

    if (species == 'Human') {
      kang = SingleR.CreateKangAnnotations(sc.data.gl)
      singler$other = kang$kang_annotation
    }

    singler$meta.data = list(project.name=project.name,orig.ident=orig.ident)

    if (reduce.file.size==T) {
      singler = remove.Unnecessary.Data.single(singler)
    }

    singler

  }
