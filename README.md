# AddVelocyto
Merge Velocyto Loom into one

#object: merged, Seurat object

#samples: sample id of each sample, should be same as velocyto output id

#cell.ids: suffix added to each cell barcode, when Seurat object merged

#input should be ordered properly

merged.velocity <- AddVelocity(files = c("sample1.loom", "sample2.loom", "sample3"), 
                                object = merged, 
                                samples = samples, 
                                cells = colnames(merged) , 
                                cell.ids = unique(substr(Cells(merged), 18, 20)))

merged.velocity@reductions <- merged@reductions

merged.velocity@meta.data <- merged@meta.data

seurat2anndata(obj = merged.velocity, outFile = "merged_velocity.h5ad")
