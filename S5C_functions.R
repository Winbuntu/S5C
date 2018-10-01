
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


find_hv_genes = function(count){
  
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

##################


GeneSetScore <- function (object, genes.list = NULL, genes.pool = NULL, n.bin = 25, 
                       seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name = "GeneSet", 
                       random.seed = 1) 
{
  
  
  object = exp.mat.norm
  genes.list = list(cc.genes$s.genes,cc.genes$g2m.genes)
  genes.pool = NULL
  n.bin = 25
  seed.use = 1
  ctrl.size = 100
  use.k = FALSE
  enrich.name = "Cluster" 
  random.seed = 1
  
  
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
                         match = rownames(x = object@data))
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





