#merged.velocity <- AddVelocity(files = c("sample1.loom", "sample2.loom", "sample3"), object = merged, samples = samples, cells = colnames(merged) , cell.ids = unique(substr(Cells(merged), 18, 22)))

AddVelocity <- function(
  files, 
  object = object,
  engine = 'hdf5r', 
  default.assay = 1,
  samples = samples,
  cell.ids = NULL,
  cells = NULL,
  slot = 'counts',
  min.cells = 0,
  min.features = 0,
  verbose = TRUE,
  ...) {
  DefaultAssay(object = object) <- "RNA"
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
      stop("cell ID required! set with cell.ids!")
    }else {
      colnames(s) <- paste(substr(colnames(s), nchar(samples[i]) + 2, nchar(samples[i]) + 17), cell.ids[i], sep = "-")
      colnames(u) <- paste(substr(colnames(u), nchar(samples[i]) + 2, nchar(samples[i]) + 17), cell.ids[i], sep = "-")
      colnames(a) <- paste(substr(colnames(a), nchar(samples[i]) + 2, nchar(samples[i]) + 17), cell.ids[i], sep = "-")
    }
    
    if (length(intersect(colnames(s), colnames(object))) == 0){
      warning("No cells in", file, " used!")
    }
    
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
  
  if (is.null(cells)){
    cells <- intersect(colnames(spliced), colnames(object))
  } else{
    cells <- intersect(intersect(colnames(spliced), colnames(object)), cells)
  }
  
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
  
  object[["RNA"]] <- CreateAssayObject(counts = object@assays$RNA@counts[rownames(spliced),])
  object[["spliced"]] <- CreateAssayObject(counts = spliced)
  object[["unspliced"]] <- CreateAssayObject(counts = unspliced)
  object[["ambiguous"]] <- CreateAssayObject(counts = ambiguous)
  
  return(object)
}
