###usage
###merged.velocity <- AddVelocyto(files = c("file1.loom", "file2.loom", "file3.loom", "file4.loom"), object = merged, samples = c("WT_c", "WT_t", "cKO_c", "cKO_t"), cell.ids = c("wc", "cc", "wp", "cp"))
###merged.velocity <- SCTransform(object = merged.velocity, assay = "spliced")
##merged.velocity <- RunVelocity(object = merged.velocity, deltaT = 1, kCells = 25, fit.quantile = 0.02)
###cell.colors <- colors[Idents(object = merged)]
###names(x = cell.colors) <- colnames(x = merged)
###show.velocity.on.embedding.cor(emb = Embeddings(object = merged.velocity, reduction = "umap"), vel = Tool(object = merged.velocity, 
#    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
#    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
#    do.par = FALSE, cell.border.alpha = 0.1)

AddVelocyto <- function(
  files, 
  object = object,
  engine = 'hdf5r', 
  default.assay = 1,
  samples = samples,
  cell.ids = NULL,
  slot = 'counts',
  min.cells = 0,
  min.features = 0,
  verbose = TRUE,
  ...) {
  DefaultAssay(object) <- "RNA"
  for (i in seq(1, length(files))) {
    file = files[i]
    if (verbose) {
      sink(file = stderr(), type = 'output')
      on.exit(expr = sink())
      x <- velocyto.R::read.loom.matrices(file = file, engine = engine)
    } else {
      invisible(x = capture.output(x <- velocyto.R::read.loom.matrices(
        file = file,
        engine = engine
      )))
    }
    
    s = x[["spliced"]][intersect(rownames(object), rownames(x[["spliced"]])), ]
    u = x[["unspliced"]][intersect(rownames(object), rownames(x[["unspliced"]])), ]
    a = x[["ambiguous"]][intersect(rownames(object), rownames(x[["ambiguous"]])), ]
    
    if (is.null(cell.ids)) {
      colnames(s) <- paste(substr(colnames(s), nchar(samples[i]) + 2, nchar(samples[i]) + 17), i, sep = "_")
      colnames(u) <- paste(substr(colnames(u), nchar(samples[i]) + 2, nchar(samples[i]) + 17), i, sep = "_")
      colnames(a) <- paste(substr(colnames(a), nchar(samples[i]) + 2, nchar(samples[i]) + 17), i, sep = "_")
    }else {
      colnames(s) <- paste(cell.ids[i], substr(colnames(s), nchar(samples[i]) + 2, nchar(samples[i]) + 17), sep = "_")
      colnames(u) <- paste(cell.ids[i], substr(colnames(u), nchar(samples[i]) + 2, nchar(samples[i]) + 17), sep = "_")
      colnames(a) <- paste(cell.ids[i], substr(colnames(a), nchar(samples[i]) + 2, nchar(samples[i]) + 17), sep = "_")
    }
    
    cells <- intersect(colnames(s), colnames(object))   
    if (i == 1){
      spliced = s
      unspliced = u
      ambiguous = a
    } else{
      spliced = cbind(spliced, s)
      unspliced = cbind(unspliced, u)
      ambiguous = cbind(ambiguous, a)
    }
  }

  cells <- intersect(colnames(spliced), colnames(object))
  if (length(cells) == 0) {
    stop("Cell number is zero!")
  }
  
  if (length(cells) != length(colnames(object))) {
    object <- subset(object, cells = cells)
  }  
  
  if (length(cells) != length(colnames(spliced))) {
    spliced <- spliced[, cells]
    unspliced <- unspliced[, cells]
    ambiguous <- ambiguous[, cells]
  }
  
  object[["spliced"]] <- CreateAssayObject(counts = spliced)
  object[["unspliced"]] <- CreateAssayObject(counts = unspliced)
  object[["ambiguous"]] <- CreateAssayObject(counts = ambiguous)
  
  return(object)
}
