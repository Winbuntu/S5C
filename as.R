files = list.files("E-GEOD-45719.processed.1/")

Deng = NA

for(f in files){
  
  temp = read.table(  paste("E-GEOD-45719.processed.1/",f,sep=""),header = F, stringsAsFactors = F)
  if(is.na(Deng)){
    Deng = temp
  }else{
    Deng = cbind(Deng, temp)
  }
  print(f)
}



deng = Deng[,c(1,seq(from = 4, to = 1902, by = 6))]

colnames(deng) = c("gene",files)

deng.emb = deng[,-c(94:114,282:318)]

colnames(deng)[c(94:106,282:318)]

cell.type = sapply(strsplit(colnames(deng.emb), split='_',fixed=TRUE),function(x) x[2])[-1]

cell.type[256:259] = "zy"

deng.emb.data = deng.emb[,-1]

####################
# deng.emb.data.rpm = sweep(x = deng.emb.data,MARGIN = 2,STATS = colSums(deng.emb.data),FUN = "/") * 1000000

cc.genes


# check if each element in the list has a minimal length "cutoff"
LengthCheck <- function(values, cutoff = 0) { 
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}



Add_info2 <- function (object, genes.list = NULL, genes.pool = NULL, n.bin = 25, 
          seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name = "Cluster", 
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


assignments = apply(X = genes.scores.use, MARGIN = 1, FUN = function(scores,first = "S", second = "G2M", null = "G1") {
  if (all(scores < 0)) {
    return(null)
  }
  else {
    return(c(first, second)[which(x = scores == max(scores))])
  }
})

############################################################################











#############################################################################



RegressOutResid <- function(
  object,
  vars.to.regress,
  genes.regress = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  display.progress = TRUE,
  do.par = FALSE,
  num.cores = 1
) {
  
  # 判断模型类型
  
  possible.models <- c("linear", "poisson", "negbinom")
  if (!model.use %in% possible.models) {
    stop(
      paste0(
        model.use,
        " is not a valid model. Please use one the following: ",
        paste0(possible.models, collapse = ", "),
        "."
      )
    )
  }
  
  
  
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- intersect(x = genes.regress, y = rownames(x = object@data))
  latent.data <- FetchData(object = object, vars.all = vars.to.regress)
  bin.size <- ifelse(test = model.use == 'negbinom', yes = 5, no = 100)
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  if (display.progress) {
    message(paste("Regressing out:", paste(vars.to.regress, collapse = ", ")))
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  }
  data.resid <- c()
  data.use <- object@data[genes.regress, , drop = FALSE];
  if (model.use != "linear") {
    use.umi <- TRUE
  }
  if (use.umi) {
    data.use <- object@raw.data[genes.regress, object@cell.names, drop = FALSE]
  }
  # input checking for parallel options
  if (do.par) {
    if (num.cores == 1) {
      num.cores <- detectCores() / 2
    } else if (num.cores > detectCores()) {
      num.cores <- detectCores() - 1
      warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
    }
  } else if (num.cores != 1) {
    num.cores <- 1
    warning("For parallel processing, please set do.par to TRUE.")
  }
  cl <- parallel::makeCluster(num.cores)#, outfile = "")
  # using doSNOW library because it supports progress bar update
  registerDoSNOW(cl)
  opts <- list()
  if (display.progress) {
    # define progress bar function
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
  }
  reg.mat.colnames <- c(colnames(x = latent.data), "GENE")
  fmla_str = paste0("GENE ", " ~ ", paste(vars.to.regress, collapse = "+"))
  if (model.use == "linear") {
    # In this code, we'll repeatedly regress different Y against the same X
    # (latent.data) in order to calculate residuals.  Rather that repeatedly
    # call lm to do this, we'll avoid recalculating the QR decomposition for the
    # latent.data matrix each time by reusing it after calculating it once and
    # storing it in a fastResiduals function.
    regression.mat <- cbind(latent.data, data.use[1,])
    colnames(regression.mat) <- reg.mat.colnames
    qr = lm(as.formula(fmla_str), data = regression.mat, qr = TRUE)$qr
    rm(regression.mat)
  }
  data.resid <- foreach(i = 1:max.bin, .combine = "c", .options.snow = opts) %dopar% {
    genes.bin.regress <- rownames(x = data.use)[bin.ind == i]
    gene.expr <- as.matrix(x = data.use[genes.bin.regress, , drop = FALSE])
    empty_char = character(length = dim(gene.expr)[1]) # Empty vector to reuse
    new.data <- sapply(
      X = genes.bin.regress,
      FUN = function(x) {
        # Fast path for std. linear models
        if(model.use=="linear") {
          resid <- qr.resid(qr, gene.expr[x,])
        } else {
          regression.mat <- cbind(latent.data, gene.expr[x,])
          colnames(x = regression.mat) <- reg.mat.colnames
          fmla = as.formula(fmla_str)
          resid <- switch(
            EXPR = model.use,
            'poisson' = residuals(
              object = glm(
                formula = fmla,
                data = regression.mat,
                family = "poisson"
              ),
              type = 'pearson'
            ),
            'negbinom' = NBResiduals(
              fmla = fmla,
              regression.mat = regression.mat,
              gene = x,
              return.mode = TRUE
            )
          )
        }
        if (!is.list(x = resid)) {
          resid <- list('resid' = resid, 'mode' = empty_char)
        }
        return(resid)
      }
    )
    new.data.resid <- new.data[seq.int(from = 1, to = length(x = new.data), by = 2)]
    new.data.resid = matrix(unlist(new.data.resid), nrow = length(new.data.resid[[1]]))
    colnames(x = new.data.resid) <- genes.bin.regress
    new.data.mode <- unlist(x = new.data[seq.int(from = 2, to = length(x = new.data), by = 2)])
    names(x = new.data.mode) <- genes.bin.regress
    new.data <- list('resid' = new.data.resid, 'mode' = new.data.mode)
    return(new.data)
  }
  if (display.progress) {
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)))
    close(pb)
  }
  stopCluster(cl)
  modes <- unlist(x = data.resid[seq.int(from = 2, to = length(x = data.resid), by = 2)])
  modes <- modes[modes == 'scale']
  names(x = modes) <- gsub(
    pattern = 'mode.',
    replacement = '',
    x = names(x = modes),
    fixed = TRUE
  )
  data.resid <- data.resid[seq.int(from = 1, to = length(x = data.resid), by = 2)]
  data.resid <- as.matrix(x = as.data.frame(x = data.resid))
  data.resid <- t(x = data.resid)
  if (length(x = modes)) {
    message(
      "The following genes failed with glm.nb, and fell back to scale(log(y+1))\n\t",
      paste(names(x = modes), collapse = ', ')
    )
  }
  rownames(x = data.resid) <- genes.regress
  suppressWarnings(expr = gc(verbose = FALSE))
  if (use.umi) {
    data.resid <- log1p(
      x = sweep(
        x = data.resid,
        MARGIN = 1,
        STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
        FUN = "-"
      )
    )
  }
  return(data.resid)
}



# Add_info <- function (object, genes.list = NULL, genes.pool = NULL, n.bin = 25, 
#                       seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name = "Cluster", 
#                       random.seed = 1) 
# {
#   set.seed(seed = random.seed)
#   
#   genes.old <- genes.list
#   
#   if (use.k) {
#     genes.list <- list()
#     for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
#       genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
#                                            i))
#     }
#     cluster.length <- length(x = genes.list)
#   }
#   else {
#     if (is.null(x = genes.list)) {
#       stop("Missing input gene list")
#     }
#     genes.list <- lapply(X = genes.list, FUN = function(x) {
#       return(intersect(x = x, y = rownames(x = object@data)))
#     })
#     cluster.length <- length(x = genes.list)
#   }
#   
#   
#   if (!all(LengthCheck(values = genes.list))) {
#     warning(paste("Could not find enough genes in the object from the following gene lists:", 
#                   paste(names(x = which(x = !LengthCheck(values = genes.list)))), 
#                   "Attempting to match case..."))
#     genes.list <- lapply(X = genes.old, FUN = CaseMatch, 
#                          match = rownames(x = object@data))
#   }
#   
#   
#   if (!all(LengthCheck(values = genes.list))) {
#     stop(paste("The following gene lists do not have enough genes present in the object:", 
#                paste(names(x = which(x = !LengthCheck(values = genes.list)))), 
#                "exiting..."))
#   }
#   
#   
#   if (is.null(x = genes.pool)) {
#     genes.pool = rownames(x = object@data)
#   }
#   
#   
#   data.avg <- apply(X = object@data[genes.pool, ], MARGIN = 1, 
#                     FUN = mean)
#   
#   
#   data.avg <- data.avg[order(data.avg)]
#   
#   
#   data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg)/n.bin)))
#   
#   
#   names(x = data.cut) <- names(x = data.avg)
#   
#   ctrl.use <- vector("list", cluster.length)
#   
#   
#   for (i in 1:cluster.length) {
#     genes.use <- genes.list[[i]]
#     for (j in 1:length(x = genes.use)) {
#       ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
#                                                                               data.cut[genes.use[j]])], size = ctrl.size, replace = FALSE)))
#     }
#   }
#   
#   ctrl.use <- lapply(X = ctrl.use, FUN = unique)
#   
#   
#   ctrl.scores <- c()
#   
#   for (i in 1:length(ctrl.use)) {
#     genes.use <- ctrl.use[[i]]
#     ctrl.scores <- rbind(ctrl.scores, apply(X = object@data[genes.use, 
#                                                             ], MARGIN = 2, FUN = mean))
#   }
#   
#   genes.scores <- c()
#   for (i in 1:cluster.length) {
#     genes.use <- genes.list[[i]]
#     if (length(genes.use) == 1) {
#       data.use <- t(as.matrix(object@data[genes.use, ]))
#     }
#     else {
#       data.use <- object@data[genes.use, ]
#     }
#     genes.scores <- rbind(genes.scores, apply(X = data.use, 
#                                               MARGIN = 2, FUN = mean))
#   }
#   
#   genes.scores.use <- genes.scores - ctrl.scores
#   rownames(x = genes.scores.use) <- paste0(enrich.name, 1:cluster.length)
#   genes.scores.use <- t(x = as.data.frame(x = genes.scores.use))
#   object <- AddMetaData(object = object, metadata = genes.scores.use, 
#                         col.name = colnames(x = genes.scores.use))
#   return(object)
# }



