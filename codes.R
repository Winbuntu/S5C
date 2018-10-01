

find_hv_genes = function(count){
  
  I = dim(count)[1]
  
  mu = apply(count, 1, function(x){  mean(x[x != log10(1.01)])   }  )
  mu[is.na(mu)] = 0
  
  sd = apply(count, 1, function(x){  sd(x[x != log10(1.01)])   }  )
  sd[is.na(sd)] = 0
  
  cv = sd/mu
  cv[is.na(cv)] = 0

  high_var_genes = which(mu >= 1 & cv >= quantile(cv, 0.25))
  
  if(length(high_var_genes) < 500){ 
    high_var_genes = 1:I} # why? Just return the whole matrix?
  
  count_hv = count[high_var_genes, ]
  return(count_hv)
  
}


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


reduce_dim <- function(count_hv){
  
  require(rsvd)
  
  pca = rsvd::rpca(t(count_hv), k = 1000, center = TRUE, scale = FALSE) 
  
  
  
  return(pca)
}





find_neighbors = function(count_hv, labeled, Kcluster = NULL, 
                          ncores, cell_labels = NULL){
  
  require(rsvd)
  require(parallel)
  require( kernlab )
  
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
    return(list(dist_list = dist_list, clust = clust))
  }
  
  if(labeled == FALSE){
    ## dimeansion reduction
    
    count_hv = deng.emb.data.norm.hv ###### test
    labeled = F
    Kcluster = 5
    J = dim(count_hv)[2]
    ncores = 1
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
    
    return(list(dist_cells = dist_cells, clust = nbs))
  }
}

nbs2 = nbs
nbs2[is.na(nbs2)]=0

plot(mat_pcs[1,], mat_pcs[2,], col = factor(nbs2))

library(cluster)

gap_stat <- clusGap(t(mat_pcs[1:2,-outliers]), FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)


first.min.cluster.num = maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")

fviz_gap_stat(gap_stat)


##########################################

# small group not regress
# no obvious group effect, do not regress

hist(count_hv[,1])



#############################################################

library(ggplot2)

deng.emb.data.norm = read_count(exp.mat)

# deng.emb.data.norm.hv = find_hv_genes(deng.emb.data.norm)

deng.emb.data.norm.hv = deng.emb.data.norm[  intersect(c(cc.genes$s.genes, cc.genes$g2m.genes) , 
                                                       rownames(deng.emb.data.norm)) ,]

deng.emb.data.norm.hv.scale = t(scale(t(deng.emb.data.norm.hv), center = T, scale = T))
#deng.emb.data.norm.hv = deng.emb.data.norm[marrow@var.genes,]

pca.res = reduce_dim(deng.emb.data.norm.hv)

library(FactoMineR)

xx.res = FactoMineR::PCA(t(deng.emb.data.norm.hv.scale),graph = F)

xx.res = FactoMineR::PCA(t(reg.ed.mat),graph = F)


plot.dat = data.frame(pc1=pca.res$x[,c(1)],
                      pc2 = pca.res$x[,c(2)],
                      type=factor(cell.type, levels = c("zy","early2cell","mid2cell","late2cell",
                                                        "4cell","8cell","16cell","earlyblast","midblast","lateblast")))

plot.dat = data.frame(pc1=pca.res$x[,c(1)],
                      pc2 = pca.res$x[,c(2)],
                      type=marrow@meta.data$Phase)

plot.dat = data.frame(pc1=xx.res$ind$coord[,1],
                      pc2 = xx.res$ind$coord[,2],
                      type=marrow@meta.data$Phase)

ggplot(plot.dat, aes(x=pc1, y=pc2, colour = type)) + geom_point()

########################


dat.mat.reg = deng.emb.data.norm.hv.scale

reg.ed.mat = matrix(NA, nrow = nrow(dat.mat.reg),
                    ncol = ncol(dat.mat.reg))


for( i in c(1:nrow(dat.mat.reg))){
  
  exp1 = as.numeric(deng.emb.data.norm.hv[i,])
  
  lm.dat = data.frame(exp1=exp1, s.score = genes.scores.use[,1], g2m.score = genes.scores.use[,2])
  
  lmEC <- lm(exp1 ~ s.score + g2m.score, data = lm.dat)
  ## Save residuals
  
  reg.ed.mat[i,] = residuals(lmEC)
}







