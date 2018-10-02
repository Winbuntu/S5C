load("cc.geges.RDS") # this is required


read_count <- function(raw_count){
  
  totalCounts_by_cell = colSums(raw_count)
  totalCounts_by_cell[totalCounts_by_cell == 0] = 1
  raw_count = sweep(raw_count, MARGIN = 2, 10^6/totalCounts_by_cell, FUN = "*")
  
  if (min(raw_count) < 0) {
    stop("smallest read count cannot be negative!")
  }
  
  count_lnorm = log10(raw_count + 1.01)
  
  return(count_lnorm)
}


find_hv_genes <- function(count){
  
  mu = apply(count, 1, function(x){  mean(x[x != log10(1.01)])   }  )
  mu[is.na(mu)] = 0
  
  sd = apply(count, 1, function(x){  sd(x[x != log10(1.01)])   }  )
  sd[is.na(sd)] = 0
  
  cv = sd/mu
  cv[is.na(cv)] = 0
  
  high_var_genes = which(mu >= 1 & cv >= quantile(cv, 0.25))
  
  if(length(high_var_genes) < 500){ 
    return(count)
  }else{
    return(count[high_var_genes, ])
    }

}

##########

# Adopted from Seurat
LengthCheck <- function(values, cutoff = 0) { 
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}

# Adopted from Seurat 
CaseMatch <- function (search, match) 
{
  search.match <- sapply(X = search, FUN = function(s) {
    return(grep(pattern = paste0("^", s, "$"), x = match, 
                ignore.case = TRUE, perl = TRUE, value = TRUE))
  })
  return(unlist(x = search.match))
}


ScaleMatrix <- function(mat){
  
  return(t(scale(t(mat), center = T, scale = T)))
  
  
}



find_neighbors <- function(count_hv, labeled=F, Kcluster = NULL, 
                          ncores=1, cell_labels = NULL){
  
  require(rsvd)
  require(parallel)
  require(kernlab)
  
  J = dim(count_hv)[2]
  
  if(labeled == TRUE){
    if(class(cell_labels) == "character"){
      labels_uniq = unique(cell_labels)
      labels_mth = 1:length(labels_uniq)
      names(labels_mth) = labels_uniq
      clust = labels_mth[cell_labels]
    }else{
      clust = cell_labels
    }
    nclust = length(unique(clust))
    print("calculating cell distances ...")
    dist_list = lapply(1:nclust, function(ll){
      cell_inds = which(clust == ll)
      count_hv_sub = count_hv[, cell_inds, drop = FALSE]
      if(length(cell_inds) < 1000){
        var_thre = 0.4
        pca = prcomp(t(count_hv_sub))
        eigs = (pca$sdev)^2
        var_cum = cumsum(eigs)/sum(eigs)
        if(max(var_cum) <= var_thre){
          npc = length(var_cum)
        }else{
          npc = which.max(var_cum > var_thre)
          if (labeled == FALSE){ npc = max(npc, Kcluster) }
        }
      }else{
        var_thre = 0.6
        pca = rpca(t(count_hv_sub), k = 1000, center = TRUE, scale = FALSE) 
        eigs = (pca$sdev)^2
        var_cum = cumsum(eigs)/sum(eigs)
        if(max(var_cum) <= var_thre){
          npc = length(var_cum)
        }else{
          npc = which.max(var_cum > var_thre)
          if (labeled == FALSE){ npc = max(npc, Kcluster) }
        }
      }
      
      if (npc < 3){ npc = 3 }
      mat_pcs = t(pca$x[, 1:npc]) 
      
      dist_cells_list = mclapply(1:length(cell_inds), function(id1){
        d = sapply(1:id1, function(id2){
          sse = sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
          sqrt(sse)
        })
        return(c(d, rep(0, length(cell_inds)-id1)))
      }, mc.cores = ncores)
      dist_cells = matrix(0, nrow = length(cell_inds), ncol = length(cell_inds))
      for(cellid in 1:length(cell_inds)){dist_cells[cellid, ] = dist_cells_list[[cellid]]}
      dist_cells = dist_cells + t(dist_cells)
      return(dist_cells)
    })
    return(list(dist_list = dist_list, clust = clust, pc_score= t(mat_pcs)))
  }
  
  if(labeled == FALSE){
    ## dimeansion reduction
    
    #count_hv = deng.emb.data.norm.hv ###### test
    #labeled = F
    #Kcluster = 5
    # J = dim(count_hv)[2]
    #ncores = 1
    # J = 4000
    
    print("dimension reduction ...")
    if(J < 5000){
      var_thre = 0.4
      pca = prcomp(t(count_hv))
      eigs = (pca$sdev)^2
      var_cum = cumsum(eigs)/sum(eigs)
      if(max(var_cum) <= var_thre){
        npc = length(var_cum)
      }else{
        npc = which.max(var_cum > var_thre)
        if (labeled == FALSE){ npc = max(npc, Kcluster) }
      }
    }else{
      var_thre = 0.6
      pca = rpca(t(count_hv), k = 1000, center = TRUE, scale = FALSE) 
      eigs = (pca$sdev)^2
      var_cum = cumsum(eigs)/sum(eigs)
      if(max(var_cum) <= var_thre){
        npc = length(var_cum)
      }else{
        npc = which.max(var_cum > var_thre)
        if (labeled == FALSE){ npc = max(npc, Kcluster) }
      }
    }
    
    if (npc < 3){ npc = 3 }
    
    mat_pcs = t(pca$x[, 1:npc]) # columns are cells
    
    
    ## detect outliers
    print("calculating cell distances ...")
    dist_cells_list = mclapply(1:J, function(id1){
      d = sapply(1:id1, function(id2){
        sse = sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
        sqrt(sse)
      })
      return(c(d, rep(0, J-id1)))
    }, mc.cores = ncores)
    dist_cells = matrix(0, nrow = J, ncol = J)
    for(cellid in 1:J){dist_cells[cellid, ] = dist_cells_list[[cellid]]}
    dist_cells = dist_cells + t(dist_cells)
    
    min_dist = sapply(1:J, function(i){
      min(dist_cells[i, -i])
    })
    
    iqr = quantile(min_dist, 0.75) - quantile(min_dist, 0.25)
    
    outliers = which(min_dist > 1.5 * iqr + quantile(min_dist, 0.75))
    
    ## clustering
    # non_out = setdiff(1:J, outliers)
    
    spec_res = specc(t(mat_pcs[, -outliers]), centers = Kcluster, kernel = "rbfdot")
    print("cluster sizes:")
    print(spec_res@size)
    nbs = rep(NA, J)
    nbs[-outliers] = spec_res
    
    return(list(dist_cells = dist_cells, clust = nbs, pc_score= t(mat_pcs) ))
  }
}


##################


GeneSetScore <- function (object, genes.list = list(cc.genes$s.genes,cc.genes$g2m.genes), genes.pool = NULL, n.bin = 25, 
                       seed.use = 1, ctrl.size = 100, enrich.name = "GeneSet", 
                       random.seed = 1) 
{
  
  
  # object = exp.mat.norm
  # genes.list = list(cc.genes$s.genes,cc.genes$g2m.genes)
  # genes.pool = NULL
  # n.bin = 25
  # seed.use = 1
  # ctrl.size = 100
  # use.k = FALSE
  # enrich.name = "Cluster" 
  # random.seed = 1
  
  
  set.seed(seed = random.seed)
  
  genes.old <- genes.list
  
  
  if (is.null(x = genes.list)) {
    stop("Missing input gene list")
  }
  
  genes.list <- lapply(X = genes.list, FUN = function(x) { # keep only those genes in the origional data frame
    return(intersect(x = x, y = rownames(x = object)))
  })
  
  
  cluster.length <- length(x = genes.list) # number of elements in the list -> 1
  
  
  
  if (!all(LengthCheck(values = genes.list))) { 
    warning(paste("Could not find enough genes in the object from the following gene lists:", # make sure there are enough genes
                  paste(names(x = which(x = !LengthCheck(values = genes.list)))), 
                  "Attempting to match case..."))
    genes.list <- lapply(X = genes.old, FUN = CaseMatch, 
                         match = rownames(x = object))
  }
  
  
  if (!all(LengthCheck(values = genes.list))) {
    stop(paste("The following gene lists do not have enough genes present in the object:", # make sure there are enough genes
               paste(names(x = which(x = !LengthCheck(values = genes.list)))), 
               "exiting..."))
  }
  
  
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(x = object) # gene pool, all genes
  }
  
  
  data.avg <- apply(X = object[genes.pool, ], MARGIN = 1, 
                    FUN = mean)
  
  
  data.avg <- data.avg[order(data.avg)] # rank by average expression
  
  
  data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg)/n.bin)))
  
  
  names(x = data.cut) <- names(x = data.avg)
  
  ctrl.use <- vector("list", cluster.length) # 做 control list
  
  
  # 从每个bin里面 sample control genes
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(x = genes.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                              data.cut[genes.use[j]])], size = ctrl.size, replace = FALSE)))
    }
  }
  
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  
  
  ctrl.scores <- c()
  
  for (i in 1:length(ctrl.use)) { # ctrl score 是control genes 的平均表达量
    genes.use <- ctrl.use[[i]]
    ctrl.scores <- rbind(ctrl.scores, apply(X = object[genes.use, 
                                                       ], MARGIN = 2, FUN = mean))
  }
  
  
  genes.scores <- c()
  
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    
    if (length(genes.use) == 1) {
      data.use <- t(as.matrix(object[genes.use, ]))
    }
    else {
      data.use <- object[genes.use, ]
    } # 确保只有一个 gene的时候取出来不会出错
    
    genes.scores <- rbind(genes.scores, apply(X = data.use, 
                                              MARGIN = 2, FUN = mean))
  }
  
  genes.scores.use <- genes.scores - ctrl.scores
  
  rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
  
  
  genes.scores.use <- t(x = as.data.frame(x = genes.scores.use))
  
  return(genes.scores.use)
  #object <- AddMetaData(object = object, metadata = genes.scores.use, 
  #                      col.name = colnames(x = genes.scores.use))
  #return(object)
}

##################
# regression

RegressOut_clusters <- function(dat.mat.reg, genes.scores.use, cluster_indicater, minimal_to_regress = 20){
  
  dat.mat.reg = data.matrix(deng.emb.dedup.data.norm)
  
  cluster_indicater = deng.emb.dedup.data.hv.clusters$clust
  minimal_to_regress = 20

  reg.ed.mat = matrix(NA, nrow = nrow(dat.mat.reg),
                      ncol = ncol(dat.mat.reg))
  
  
  
  for(l in  na.omit(unique(cluster_indicater))   ){
    
    print(l)
    
    # l = 4
    # skip regression if have too few cells
    
    print(dim(reg.ed.mat))
    
  if(length(which(cluster_indicater==l)) < minimal_to_regress ){
    
    reg.ed.mat[, which(    cluster_indicater==l   )] = dat.mat.reg[, which(  cluster_indicater==l   )]
    
    print(paste0("skip cluster ", l))
    next
  }
  
    
    this.level.exp = dat.mat.reg[,  which(cluster_indicater==l)  ]
    this.level.genes.scores.use = GeneSetScore(object = this.level.exp)
    
    
  for( i in c(1:nrow(dat.mat.reg))){
    
    #i=1
    exp1 = as.numeric(dat.mat.reg[i,  which(cluster_indicater==l)  ])
    
    
    
    lm.dat = data.frame(exp1=exp1, s.score = this.level.genes.scores.use[ ,1], 
                        g2m.score = this.level.genes.scores.use[ ,2])
    
    lmEC <- lm(exp1 ~ s.score + g2m.score, data = lm.dat)
    ## Save residuals
    
    reg.ed.mat[i, which(cluster_indicater==l)] = residuals(lmEC)

  }
  
  }
  
  # fill outliers
  
  reg.ed.mat[,which(is.na( cluster_indicater ))] = dat.mat.reg[,which(is.na( cluster_indicater ))]
  
  return(reg.ed.mat)
}


#############################################

head(deng.emb)

deng.emb.dedup = deng.emb[(!duplicated(deng.emb$gene)),]
deng.emb.dedup.data = deng.emb.dedup
rownames(deng.emb.dedup.data) = deng.emb.dedup$gene
deng.emb.dedup.data = deng.emb.dedup.data[,-1]

################

deng.emb.dedup.data.norm = read_count(deng.emb.dedup.data)

deng.emb.dedup.data.hv = find_hv_genes(deng.emb.dedup.data.norm)

deng.emb.dedup.data.hv.scaled = ScaleMatrix(deng.emb.dedup.data.hv)

deng.emb.dedup.data.hv.scaled.pca.res = FactoMineR::PCA(t(deng.emb.dedup.data.hv),graph = F)


plot.dat = data.frame(pc1=deng.emb.dedup.data.hv.scaled.pca.res$ind$coord[,c(1)],
                      pc2 = deng.emb.dedup.data.hv.scaled.pca.res$ind$coord[,c(2)],
                      type=factor(cell.type, levels = c("zy","early2cell","mid2cell","late2cell",
                                                        "4cell","8cell","16cell","earlyblast","midblast","lateblast")))
ggplot(plot.dat, aes(x=pc1, y=pc2, colour = type)) + geom_point()


########## find clusters ########

deng.emb.dedup.data.hv.clusters = find_neighbors(count_hv = deng.emb.dedup.data.hv, labeled = F,Kcluster = 5,ncores = 1)

plot.dat.spec.cluster = data.frame(pc1=deng.emb.dedup.data.hv.scaled.pca.res$ind$coord[,c(1)],
                      pc2 = deng.emb.dedup.data.hv.scaled.pca.res$ind$coord[,c(2)],
                      type=factor(deng.emb.dedup.data.hv.clusters$clust)  )

ggplot(plot.dat.spec.cluster, aes(x=pc1, y=pc2, colour = type)) + geom_point()


##############################

# 对每一个cluster regress out



deng.emb.dedup.data.ccc = RegressOut_clusters(dat.mat.reg = deng.emb.dedup.data.norm, cluster_indicater = deng.emb.dedup.data.hv.clusters$clust )











