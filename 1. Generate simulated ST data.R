########simulated st data based on sc data#########
######use mean as the value for simulated spot#######

library(rmutil)
library(gtools)
library(Matrix)
library(Seurat)
library(dplyr)
library(data.table)
#########################
scmeta<-read.table("./simulation/GSE182227_OPSCC/metadata.txt",header=T,sep="\t")
counts <- as.matrix(readMM("./simulation/GSE182227_OPSCC/matrix.mtx.gz"))
genes <- read.table("./simulation/GSE182227_OPSCC/features.tsv", sep = '\t', header = FALSE)
barcodes <- read.table("./simulation/GSE182227_OPSCC/barcodes.tsv", sep = '\t', header = FALSE)

# Annotate rows and columns of counts matrix
rownames(counts) <- genes[, 1]
colnames(counts) <- barcodes[, 1]

# Create a Seurat object
seurat_object <- CreateSeuratObject(counts = counts)
n_spots <- 1000
min_cells_per_spot <- 3
max_cells_per_spot <- 10
max_cell_types_per_spot <- 3
for (sam in c("OP5","OP8","OP10","OP12","OP13")){
print(sam)
# M_spots_count <- round(M_appearance_rates[[sam]] * n_spots)
sc_object <-subset(seurat_object,subset =orig.ident== sam)

n_cells <- 10
test_scmeta<-scmeta%>%filter(Patient==sam)
test_scmeta<-test_scmeta%>%filter(Type!="Unresolved")
test_scmeta$Type<-ifelse(test_scmeta$Cancer_assignment=="Cancer","Malignant",test_scmeta$Type)
rm_type<-as.data.frame(table(test_scmeta$Type))%>%filter(Freq<10)%>%.$Var1%>%as.character()
test_scmeta<-test_scmeta%>%filter(!Type %in% rm_type)

#test_scmeta<-test_scmeta[which(test_scmeta$nCount_RNA>1000),]
cell_types <-unique(test_scmeta$Type)
cell_data <-sc_object$RNA@data


# Create an empty list to store the simulated spots
spots <- vector("list", n_spots)
spots_celltype_info <- vector("list", n_spots)
# Simulate each spot
cell_type_table <- matrix(NA, nrow = n_spots, ncol = 6, 
                          dimnames = list(NULL, c("CellType1", "CellType2", "CellType3", "Prop1", "Prop2", "Prop3")))
set.seed(123)  
for (i in 1:n_spots) {
 # Randomly decide the number of cells in this spot
  n_cells_in_spot <- sample(min_cells_per_spot:max_cells_per_spot, 1)  
  # Randomly decide the number of cell types in this spot, not exceeding the max_cell_types_per_spot
  n_cell_types_in_spot <- sample(1:max_cell_types_per_spot, 1)
  #print(n_cell_types_in_spot)
  # Randomly decide the cell types in this spot
  
#########
  cell_types_in_spot <- sample(cell_types, n_cell_types_in_spot, replace = FALSE)
  #print(paste0("cell_type:",length(unique(cell_types_in_spot))))
  # Generate cell type proportions using the Dirichlet distribution
  cell_type_proportions <- rdirichlet(1, rep(1, length(cell_types_in_spot)))[1,] 
  cell_type_table[i, 1:length(cell_types_in_spot)] <- cell_types_in_spot
  cell_type_table[i, (1:length(cell_type_proportions)) + 3] <- cell_type_proportions
  # For each cell type, sample the corresponding number of cells from the original data
  spot_cells <- lapply(1:n_cells_in_spot, function(x) {
    cell_type <- sample(cell_types_in_spot, 1, prob = cell_type_proportions)
    test_scmeta[sample(which(test_scmeta$Type == cell_type), 1), ]
  })
  spots_celltype_info[[i]] <- sapply(spot_cells, function(cell) cell$Type)
  # Combine the cells into a spot using the mean
  spots[[i]] <- apply(cell_data[,sapply(spot_cells, function(cell) cell$Cell_barcode)],1,mean)
   # Combine the cells into a spot using the sum
  #spots[[i]] <- apply(cell_data[,sapply(spot_cells, function(cell) cell$Cell_barcode)],1,sum)
}
rownames(cell_type_table)<-paste0("spot", 1:nrow(cell_type_table))
write.table(cell_type_table,paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type_propotion_sum.txt"),sep="\t")
  
simulated_data <- do.call(rbind, spots)
rownames(simulated_data) <- paste0("spot", 1:nrow(simulated_data))
simulated_data<-t(simulated_data)
# Combine the cell type info into a data frame
spots_celltype_info <- lapply(spots_celltype_info, function(x) paste(x, collapse = ","))
# Convert the list into a dataframe
spots_celltype_info_df <- data.frame(Spot = 1:n_spots,CellTypes = unlist(spots_celltype_info))
spots_celltype_info_df$Num_Cell_Types <- sapply(spots_celltype_info, function(x) length(unique(strsplit(x, ",")[[1]])))
spots_celltype_info_df$Num_Cells <- sapply(spots_celltype_info, function(x) length(strsplit(x, ",")[[1]]))
rownames(spots_celltype_info_df) <- paste0("spot", 1:nrow(spots_celltype_info_df))
print(table(spots_celltype_info_df$Num_Cell_Types))
write.table(simulated_data,paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/",sam,"si_ST_mtx_sum.txt"),sep="\t")
write.table(spots_celltype_info_df,paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type_sum.txt"),sep="\t")
}




simulated_data <- simulated_data[rowSums(simulated_data) > 0, ]

# Filter columns where sum is greater than zero
simulated_data <- simulated_data[ , colSums(simulated_data) > 0]

range(apply(simulated_data,2,function(x)length(which(x>0))))
range(rowSums(simulated_data))
range(colSums(simulated_data))
gsva_res<-data.frame()
for (sam in c("OP5","OP8","OP10","OP12","OP13")){
print(sam)
SCTdata<-read.table(paste0("./simulation/GSE182227_OPSCC/si_ST/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
SCTdata<-as.matrix(SCTdata)
print(range(apply(SCTdata,2,function(x)length(which(x>0)))))
print(range(rowSums(SCTdata)))
print(range(colSums(SCTdata)))


gsva_score<-as.data.frame(t(gsva(SCTdata, markerlist, parallel.sz=1)))
gsva_res<-rbind(gsva_res,gsva_score)
write.table(gsva_res,"./simulation/GSE182227_OPSCC/si_ST/gsva_htca_cancer_res_6.txt",sep="\t")
}

get_frequencies <- function(vector_string) {
  vector_elements <- strsplit(as.character(vector_string), split = ",")[[1]]
  return(table(vector_elements))
}

for (sam in c("OP5","OP8","OP10","OP12","OP13")){
print(sam)
spotdata<-read.table(paste0("./simulation/GSE182227_OPSCC/si_ST/",sam,"si_ST_cell_type.txt"),sep="\t",header=T)
SCTdata<-read.table(paste0("./simulation/GSE182227_OPSCC/si_ST/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
frequencies <- lapply(spotdata$CellTypes, get_frequencies)
# Initialize a list to store the frequency data for each row
frequency_list <- list()
# Iterate through the frequencies and convert to data frame
for (i in seq_along(frequencies)) {
  frequency_list[[i]] <- as.data.frame(as.table(frequencies[[i]]))
}
frequency_df <- data.frame(element = character(), frequency = integer(),  percent=numeric(),row_name = character())
# Iterate through the frequency list, adding the required information to the final data frame
for (i in seq_along(frequency_list)){
  freq_table <- frequency_list[[i]]
  elements <- freq_table$vector_elements
  frequencies <- freq_table$Freq
  percentages <- frequencies / sum(frequencies) 
  row_df <- data.frame(
    element = elements,
    frequency = frequencies,
	percent=percentages,
    row_name = rep(rownames(spotdata)[i], length(elements)), # Repeat the row name for each element
    stringsAsFactors = FALSE
  )
  frequency_df <- rbind(frequency_df, row_df)
}
frequency_df$dataset<-sam
nor_spot<- frequency_df %>%
  group_by(row_name) %>%
  filter(!any(element == "Malignant")) %>%
  select(row_name) %>%
  distinct()
set.seed(222)
nor_spots <- sample(nor_spot$row_name, 100, replace = FALSE)
nor_spots<-paste0(sam,"_",nor_spots)
nor_data<-SCTdata[,nor_spots]
write.table(frequency_df,paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type_summmary.txt"),sep="\t")
write.table(nor_data,paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/",sam,"si_ST_normal_matrix.txt"),sep="\t")
}


nor_matrix<-read.table(paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/","OP5","si_ST_normal_matrix.txt"),sep="\t")
nor_matrix$gene<-rownames(nor_matrix)
for (sam in c("OP8","OP10","OP12","OP13")){
print(sam)
nor_data<-read.table(paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/",sam,"si_ST_normal_matrix.txt"),sep="\t")
nor_data$gene<-rownames(nor_data)
nor_matrix<-merge(nor_matrix,nor_data,by="gene")
}
rownames(nor_matrix)<-nor_matrix$gene
nor_matrix<-nor_matrix[,-1]

write.table(nor_matrix,"./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/si_ST_normal_matrix.txt",sep="\t")

normal_count<-nor_matrix
cell.data<-as.data.frame(colnames(normal_count))
cell.data$sample_info<-sapply(strsplit(cell.data[,1], "_"), function(x) paste(x[-length(x)], collapse = "_"))
colnames(cell.data)[1]<-"cell.names"


for (sam in c("OP5","OP8","OP10","OP12","OP13")){
print(sam)
SCTdata<-read.table(paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
norm.cell.names<-cell.data$cell.names[which(cell.data$sample_info==sam)]
spatial_count<-SCTdata[,which(!(colnames(SCTdata) %in% norm.cell.names))]
gene_to_keep <- apply(spatial_count, 1, function(row) length(which(row > 0)) >= 100)
spatial_filtered <-spatial_count[gene_to_keep,]
rows_to_remove <- grepl("^MT-", rownames(spatial_filtered))
# Subset the dataframe to keep only the rows that don't start with "MT-"
spatial_filtered <- spatial_filtered[!rows_to_remove,]
print(dim(spatial_filtered))
spatial_filtered$gene<-rownames(spatial_filtered) 
normal_count$gene<-rownames(normal_count) 
all_count<-merge(spatial_filtered,normal_count,by="gene")
rownames(all_count)<-all_count$gene
all_count<-all_count[,-1]
print(length(intersect(colnames(all_count),cell.data$cell.names)))

print("start cnv")
cnv_data<-copykat(rawmat=all_count,id.type="S", LOW.DR=0.03, UP.DR=0.2,ngene.chr=5,sam.name=paste0(sam,"_st_1017_pan"),KS.cut=0.1,win.size=25, norm.cell.names= cell.data$cell.names)
}

cnv_res<-data.frame()
for (sam in c("OP5","OP8","OP10","OP12","OP13")){
print(sam)

cnv_data<-read.table(paste0("./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/cnv/",sam,"_st_1017_copykat_prediction.txt"),header=T,sep="\t")
cnv_data<-cnv_data%>%filter(copykat.pred!="not.defined"& grepl(paste0("^",sam), cell.names))
cnv_res<-rbind(cnv_res,cnv_data)
}
write.table(cnv_res,"./simulation/GSE182227_OPSCC/si_ST/dif_Mspots_percet/cnv/5_sample_cnv_res.txt",sep="\t")


###################
scmeta<-read.table("./simulation/GSE165897_OVCA/GSE165897_cellInfo_HGSOC.tsv",header=T,sep="\t")
counts <- fread("./simulation/GSE165897_OVCA/GSE165897_UMIcounts_HGSOC.tsv.gz")
counts2<-as.matrix(counts[,2:ncol(counts)])
rownames(counts2)<-as.character(counts$V1)

seurat_object <- CreateSeuratObject(counts = counts2)
seurat_object$orig.ident<-sapply(strsplit(rownames(seurat_object@meta.data), "-"),function(x) unlist(x)[2])

n_spots <- 1000
min_cells_per_spot <- 3
max_cells_per_spot <- 10
max_cell_types_per_spot <- 3


for (sam in  1:length(unique(scmeta$sample))){
print(sam)
#M_spots_count <- round(M_appearance_rates[[sam]] * n_spots)
sc_object <-subset(seurat_object,subset =orig.ident== unique(seurat_object$orig.ident)[sam])
print(paste0(unique(seurat_object$orig.ident)[sam],",",unique(scmeta$sample)[sam]))
n_cells <- 10
test_scmeta<-scmeta%>%filter(sample==unique(scmeta$sample)[sam])
#test_scmeta<-test_scmeta%>%filter(Type!="Unresolved")
test_scmeta$Type<-test_scmeta$cell_subtype 
test_scmeta$Type[test_scmeta$cell_subtype %in% paste0("EOC_C",1:12)] <- "Malignant"
test_scmeta$Type[test_scmeta$cell_subtype %in% paste0("CAF-",1:3)] <- "CAF"
test_scmeta$Type[test_scmeta$cell_subtype %in% paste0("DC-",1:2)] <- "DC"

rm_type<-as.data.frame(table(test_scmeta$Type))%>%filter(Freq<10)%>%.$Var1%>%as.character()
test_scmeta<-test_scmeta%>%filter(!Type %in% rm_type)

#test_scmeta<-test_scmeta[which(test_scmeta$nCount_RNA>1000),]
cell_types <-unique(test_scmeta$Type)
cell_data <-sc_object$RNA@data

# Create an empty list to store the simulated spots
spots <- vector("list", n_spots)
spots_celltype_info <- vector("list", n_spots)
# Simulate each spot
cell_type_table <- matrix(NA, nrow = n_spots, ncol = 6, 
                          dimnames = list(NULL, c("CellType1", "CellType2", "CellType3", "Prop1", "Prop2", "Prop3")))
set.seed(123)  
for (i in 1:n_spots) {
 # Randomly decide the number of cells in this spot
  n_cells_in_spot <- sample(min_cells_per_spot:max_cells_per_spot, 1)  
  # Randomly decide the number of cell types in this spot, not exceeding the max_cell_types_per_spot
  n_cell_types_in_spot <- sample(1:max_cell_types_per_spot, 1)
  #print(n_cell_types_in_spot)
  # Randomly decide the cell types in this spot
	
  cell_types_in_spot <- sample(cell_types, n_cell_types_in_spot, replace = FALSE)
  #print(paste0("cell_type:",length(unique(cell_types_in_spot))))
  # Generate cell type proportions using the Dirichlet distribution
  cell_type_proportions <- rdirichlet(1, rep(1, length(cell_types_in_spot)))[1,] 
  cell_type_table[i, 1:length(cell_types_in_spot)] <- cell_types_in_spot
  cell_type_table[i, (1:length(cell_type_proportions)) + 3] <- cell_type_proportions
  # For each cell type, sample the corresponding number of cells from the original data
  spot_cells <- lapply(1:n_cells_in_spot, function(x) {
    cell_type <- sample(cell_types_in_spot, 1, prob = cell_type_proportions)
    test_scmeta[sample(which(test_scmeta$Type == cell_type), 1), ]
  })
  spots_celltype_info[[i]] <- sapply(spot_cells, function(cell) cell$Type)
  # Combine the cells into a spot using the mean
  spots[[i]] <- apply(cell_data[,sapply(spot_cells, function(cell) cell$cell)],1,mean)
}
rownames(cell_type_table)<-paste0("spot", 1:nrow(cell_type_table))
write.table(cell_type_table,paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/", unique(seurat_object$orig.ident)[sam],"si_ST_cell_type_propotion.txt"),sep="\t")
  
simulated_data <- do.call(rbind, spots)
rownames(simulated_data) <- paste0("spot", 1:nrow(simulated_data))
simulated_data<-t(simulated_data)
# Combine the cell type info into a data frame
spots_celltype_info <- lapply(spots_celltype_info, function(x) paste(x, collapse = ","))
# Convert the list into a dataframe
spots_celltype_info_df <- data.frame(Spot = 1:n_spots,CellTypes = unlist(spots_celltype_info))
spots_celltype_info_df$Num_Cell_Types <- sapply(spots_celltype_info, function(x) length(unique(strsplit(x, ",")[[1]])))
spots_celltype_info_df$Num_Cells <- sapply(spots_celltype_info, function(x) length(strsplit(x, ",")[[1]]))
rownames(spots_celltype_info_df) <- paste0("spot", 1:nrow(spots_celltype_info_df))
print(table(spots_celltype_info_df$Num_Cell_Types))
write.table(simulated_data,paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/", unique(seurat_object$orig.ident)[sam],"si_ST_mtx.txt"),sep="\t")
write.table(spots_celltype_info_df,paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/", unique(seurat_object$orig.ident)[sam],"si_ST_cell_type.txt"),sep="\t")
}
simulated_data <- simulated_data[rowSums(simulated_data) > 0, ]

# Filter columns where sum is greater than zero
simulated_data <- simulated_data[ , colSums(simulated_data) > 0]

range(apply(simulated_data,2,function(x)length(which(x>0))))
range(rowSums(simulated_data))
range(colSums(simulated_data))

gsva_res<-data.frame()
for (sam in  1:length(unique(scmeta$sample))){
print(paste0(unique(seurat_object$orig.ident)[sam],",",unique(scmeta$sample)[sam]))
SCTdata<-read.table(paste0("./simulation/GSE165897_OVCA/si_ST/",unique(seurat_object$orig.ident)[sam],"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(unique(seurat_object$orig.ident)[sam],"_",colnames(SCTdata))
SCTdata<-as.matrix(SCTdata)
print(range(apply(SCTdata,2,function(x)length(which(x>0)))))
print(range(rowSums(SCTdata)))
print(range(colSums(SCTdata)))
gsva_score<-as.data.frame(t(gsva(SCTdata, markerlist, parallel.sz=1)))
gsva_res<-rbind(gsva_res,gsva_score)
write.table(gsva_res,"./simulation/GSE165897_OVCA/si_ST/gsva_htca_cancer_res_6.txt",sep="\t")
}

get_frequencies <- function(vector_string) {
  vector_elements <- strsplit(as.character(vector_string), split = ",")[[1]]
  return(table(vector_elements))
}

for (sam in 1:length(unique(scmeta$sample))){
print(paste0(unique(seurat_object$orig.ident)[sam],",",unique(scmeta$sample)[sam]))
spotdata<-read.table(paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/",unique(seurat_object$orig.ident)[sam],"si_ST_cell_type.txt"),sep="\t",header=T)
SCTdata<-read.table(paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/",unique(seurat_object$orig.ident)[sam],"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(unique(seurat_object$orig.ident)[sam],"_",colnames(SCTdata))
frequencies <- lapply(spotdata$CellTypes, get_frequencies)
# Initialize a list to store the frequency data for each row
frequency_list <- list()
# Iterate through the frequencies and convert to data frame
for (i in seq_along(frequencies)) {
  frequency_list[[i]] <- as.data.frame(as.table(frequencies[[i]]))
}
frequency_df <- data.frame(element = character(), frequency = integer(),  percent=numeric(),row_name = character())
# Iterate through the frequency list, adding the required information to the final data frame
for (i in seq_along(frequency_list)){
  freq_table <- frequency_list[[i]]
  elements <- freq_table$vector_elements
  frequencies <- freq_table$Freq
  percentages <- frequencies / sum(frequencies) 
  row_df <- data.frame(
    element = elements,
    frequency = frequencies,
	percent=percentages,
    row_name = rep(rownames(spotdata)[i], length(elements)), # Repeat the row name for each element
    stringsAsFactors = FALSE
  )
  frequency_df <- rbind(frequency_df, row_df)
}
frequency_df$dataset<-unique(seurat_object$orig.ident)[sam]
nor_spot<- frequency_df %>%
  group_by(row_name) %>%
  filter(!any(element == "Malignant")) %>%
  select(row_name) %>%
  distinct()
set.seed(222)
nor_spots <- sample(nor_spot$row_name, 10, replace = FALSE)
nor_spots<-paste0(unique(seurat_object$orig.ident)[sam],"_",nor_spots)
nor_data<-SCTdata[,nor_spots]
write.table(frequency_df,paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/",unique(seurat_object$orig.ident)[sam],"si_ST_cell_type_summmary.txt"),sep="\t")
write.table(nor_data,paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/",unique(seurat_object$orig.ident)[sam],"si_ST_normal_matrix.txt"),sep="\t")
}


nor_matrix<-read.table(paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/",unique(seurat_object$orig.ident)[1],"si_ST_normal_matrix.txt"),sep="\t")
nor_matrix$gene<-rownames(nor_matrix)
for (sam in 2:length(unique(scmeta$sample))){
print(paste0(unique(seurat_object$orig.ident)[sam],",",unique(scmeta$sample)[sam]))
nor_data<-read.table(paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/",unique(seurat_object$orig.ident)[sam],"si_ST_normal_matrix.txt"),sep="\t")
nor_data$gene<-rownames(nor_data)
nor_matrix<-merge(nor_matrix,nor_data,by="gene")
}
rownames(nor_matrix)<-nor_matrix$gene
nor_matrix<-nor_matrix[,-1]

write.table(nor_matrix,"./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/si_ST_normal_matrix.txt",sep="\t")

normal_count<-nor_matrix
cell.data<-as.data.frame(colnames(normal_count))
cell.data$sample_info<-sapply(strsplit(cell.data[,1], "_"), function(x) paste(x[-length(x)], collapse = "_"))
colnames(cell.data)[1]<-"cell.names"


for (sam in 1:length(unique(scmeta$sample))){
print(paste0(unique(seurat_object$orig.ident)[sam],",",unique(scmeta$sample)[sam]))
SCTdata<-read.table(paste0("./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/",unique(seurat_object$orig.ident)[sam],"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(unique(seurat_object$orig.ident)[sam],"_",colnames(SCTdata))
norm.cell.names<-cell.data$cell.names[which(cell.data$sample_info==unique(seurat_object$orig.ident)[sam])]
spatial_count<-SCTdata[,which(!(colnames(SCTdata) %in% norm.cell.names))]
gene_to_keep <- apply(spatial_count, 1, function(row) length(which(row > 0)) >= 100)
spatial_filtered <-spatial_count[gene_to_keep,]
rows_to_remove <- grepl("^MT-", rownames(spatial_filtered))
# Subset the dataframe to keep only the rows that don't start with "MT-"
spatial_filtered <- spatial_filtered[!rows_to_remove,]
print(dim(spatial_filtered))
spatial_filtered$gene<-rownames(spatial_filtered) 
normal_count$gene<-rownames(normal_count) 
all_count<-merge(spatial_filtered,normal_count,by="gene")
rownames(all_count)<-all_count$gene
all_count<-all_count[,-1]
print(length(intersect(colnames(all_count),cell.data$cell.names)))

print("start cnv")
cnv_data<-copykat(rawmat=all_count,id.type="S", LOW.DR=0.03, UP.DR=0.2,ngene.chr=5,sam.name=paste0(unique(seurat_object$orig.ident)[sam],"_st_1017_pan"),KS.cut=0.1,win.size=25, norm.cell.names= cell.data$cell.names)
}

cnv_res<-data.frame()
for (sam in 1:length(unique(scmeta$sample))){
print(sam)
print(paste0(unique(seurat_object$orig.ident)[sam],",",unique(scmeta$sample)[sam]))
cnv_data<-read.table(paste0("./simulation/GSE165897_OVCA/cnv/",unique(seurat_object$orig.ident)[sam],"_st_1017_copykat_prediction.txt"),header=T,sep="\t")
cnv_data<-cnv_data%>%filter(copykat.pred!="not.defined"& grepl(paste0("^",unique(seurat_object$orig.ident)[sam]), cell.names))
cnv_res<-rbind(cnv_res,cnv_data)
}
write.table(cnv_res,"./simulation/GSE165897_OVCA/si_ST/dif_Mspots_percet/cnv/5_sample_cnv_res.txt",sep="\t")

############################################
scmeta<-read.table("./simulation/GSE224630_renal/GSE224630_overall_metadata.tsv",header=T,sep="\t")
scmeta<-scmeta%>%filter(all.cell.association!="Tubular epithelial & malignant cells")
scmeta$all.cell.association <-gsub(" ","_",scmeta$all.cell.association)

n_spots <- 1000
min_cells_per_spot <- 3
max_cells_per_spot <- 10
max_cell_types_per_spot <- 3
for (sam in c("RCC1","RCC2","RCC3","RCC4","RCC4f","RCC5")){
print(sam)
M_spots_count <- round(M_appearance_rates[[sam]] * n_spots)

seurat_object<-Read10X(paste0("./simulation/GSE224630_renal/GSE224630_RAW/",sam,"/"))
seurat_object<-CreateSeuratObject(seurat_object,min.cells = 3)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 100  &  nCount_RNA>1000& percent.mt < 20)
n_cells <- 10
cell_data <-seurat_object $RNA@data
test_scmeta<-scmeta%>%filter(sample.ID==sam)
rm_type<-as.data.frame(table(test_scmeta$all.cell.association))%>%filter(Freq<10)%>%.$Var1%>%as.character()
test_scmeta<-test_scmeta%>%filter(!all.cell.association %in% rm_type)
last_letter<-substr(test_scmeta$barcode[1], nchar(test_scmeta$barcode[1]), nchar(test_scmeta$barcode[1]))
test_scmeta$barcode<-gsub(last_letter,"1",test_scmeta$barcode)
test_scmeta<-test_scmeta[which(test_scmeta$barcode%in% colnames(cell_data)),]
#test_scmeta<-test_scmeta[which(test_scmeta$nCount_RNA>1000),]
cell_types <-unique(test_scmeta$all.cell.association)

# Create an empty list to store the simulated spots
spots <- vector("list", n_spots)
spots_celltype_info <- vector("list", n_spots)
# Simulate each spot
cell_type_table <- matrix(NA, nrow = n_spots, ncol = 6, 
                          dimnames = list(NULL, c("CellType1", "CellType2", "CellType3", "Prop1", "Prop2", "Prop3")))
set.seed(123)  
for (i in 1:n_spots) {
  #print(paste("Iteration:", i))
 # Randomly decide the number of cells in this spot
  n_cells_in_spot <- sample(min_cells_per_spot:max_cells_per_spot, 1)  
  # Randomly decide the number of cell types in this spot, not exceeding the max_cell_types_per_spot
  n_cell_types_in_spot <- sample(1:max_cell_types_per_spot, 1)
  #print(paste("n_cell_types_in_spot:", n_cell_types_in_spot))
  # Randomly decide the cell types in this spot
   
  cell_types_in_spot <- sample(cell_types, n_cell_types_in_spot, replace = FALSE)
  #print(paste0("cell_type:",length(unique(cell_types_in_spot))))
  # Generate cell type proportions using the Dirichlet distribution
  cell_type_proportions <- rdirichlet(1, rep(1, length(cell_types_in_spot)))[1,] 
  cell_type_table[i, 1:length(cell_types_in_spot)] <- cell_types_in_spot
  cell_type_table[i, (1:length(cell_type_proportions)) + 3] <- cell_type_proportions
  # For each cell type, sample the corresponding number of cells from the original data
  spot_cells <- lapply(1:n_cells_in_spot, function(x) {
    cell_type <- sample(cell_types_in_spot, 1, prob = cell_type_proportions)
    test_scmeta[sample(which(test_scmeta$all.cell.association == cell_type), 1), ]
  })
  valid_barcodes <- colnames(cell_data)
  spot_cells2 <- lapply(spot_cells, function(cell_df) {
  return(cell_df[cell_df$barcode %in% valid_barcodes, ])})

# Remove empty list elements (if any)
  spot_cells2 <- spot_cells2[sapply(spot_cells2, function(cell_df) nrow(cell_df) > 0)]
  spots_celltype_info[[i]] <- sapply(spot_cells2, function(cell) cell$all.cell.association)
  #print(spots_celltype_info[[i]])
  # Combine the cells into a spot using the mean
  
  spots[[i]] <- apply(cell_data[,sapply(spot_cells2, function(cell) cell$barcode)],1,mean)
 #print(spots[[i]][1:10])
}
rownames(cell_type_table)<-paste0("spot", 1:nrow(cell_type_table))
write.table(cell_type_table,paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type_propotion.txt"),sep="\t")
  
simulated_data <- do.call(rbind, spots)
rownames(simulated_data) <- paste0("spot", 1:nrow(simulated_data))
simulated_data<-t(simulated_data)
# Combine the cell type info into a data frame
spots_celltype_info <- lapply(spots_celltype_info, function(x) paste(x, collapse = ","))
# Convert the list into a dataframe
spots_celltype_info_df <- data.frame(Spot = 1:n_spots,CellTypes = unlist(spots_celltype_info))
spots_celltype_info_df$Num_Cell_Types <- sapply(spots_celltype_info, function(x) length(unique(strsplit(x, ",")[[1]])))
spots_celltype_info_df$Num_Cells <- sapply(spots_celltype_info, function(x) length(strsplit(x, ",")[[1]]))
rownames(spots_celltype_info_df) <- paste0("spot", 1:nrow(spots_celltype_info_df))
print(table(spots_celltype_info_df$Num_Cell_Types))
write.table(simulated_data,paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/",sam,"si_ST_mtx.txt"),sep="\t")
write.table(spots_celltype_info_df,paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type.txt"),sep="\t")
}

gsva_res<-data.frame()
for (sam in c("RCC1","RCC2","RCC3","RCC4","RCC4f","RCC5")){
print(sam)
SCTdata<-read.table(paste0("./simulation/GSE224630_renal/si_ST/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
SCTdata<-as.matrix(SCTdata)
print(range(apply(SCTdata,2,function(x)length(which(x>0)))))
print(range(rowSums(SCTdata)))
print(range(colSums(SCTdata)))
gsva_score<-as.data.frame(t(gsva(SCTdata, markerlist, parallel.sz=1)))
gsva_res<-rbind(gsva_res,gsva_score)
write.table(gsva_res,"./simulation/GSE224630_renal/si_ST/gsva_htca_cancer_res_6.txt",sep="\t")
}

get_frequencies <- function(vector_string) {
  vector_elements <- strsplit(as.character(vector_string), split = ",")[[1]]
  return(table(vector_elements))
}

for (sam in c("RCC1","RCC2","RCC3","RCC4","RCC4f","RCC5")){
print(sam)
spotdata<-read.table(paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type.txt"),sep="\t",header=T)
SCTdata<-read.table(paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
frequencies <- lapply(spotdata$CellTypes, get_frequencies)
# Initialize a list to store the frequency data for each row
frequency_list <- list()
# Iterate through the frequencies and convert to data frame
for (i in seq_along(frequencies)) {
  frequency_list[[i]] <- as.data.frame(as.table(frequencies[[i]]))
}
frequency_df <- data.frame(element = character(), frequency = integer(),  percent=numeric(),row_name = character())
# Iterate through the frequency list, adding the required information to the final data frame
for (i in seq_along(frequency_list)){
  freq_table <- frequency_list[[i]]
  elements <- freq_table$vector_elements
  frequencies <- freq_table$Freq
  percentages <- frequencies / sum(frequencies) 
  row_df <- data.frame(
    element = elements,
    frequency = frequencies,
	percent=percentages,
    row_name = rep(rownames(spotdata)[i], length(elements)), # Repeat the row name for each element
    stringsAsFactors = FALSE
  )
  frequency_df <- rbind(frequency_df, row_df)
}
frequency_df$dataset<-sam
nor_spot<- frequency_df %>%
  group_by(row_name) %>%
  filter(!any(element == "Malignant")) %>%
  select(row_name) %>%
  distinct()
set.seed(222)
nor_spots <- sample(nor_spot$row_name, 20, replace = FALSE)
nor_spots<-paste0(sam,"_",nor_spots)
nor_data<-SCTdata[,nor_spots]
write.table(frequency_df,paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type_summmary.txt"),sep="\t")
write.table(nor_data,paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/",sam,"si_ST_normal_matrix.txt"),sep="\t")
}
nor_matrix<-read.table(paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/","RCC1","si_ST_normal_matrix.txt"),sep="\t")
nor_matrix$gene<-rownames(nor_matrix)
for (sam in c("RCC2","RCC3","RCC4","RCC4f","RCC5")){
print(sam)
nor_data<-read.table(paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/",sam,"si_ST_normal_matrix.txt"),sep="\t")
nor_data$gene<-rownames(nor_data)
nor_matrix<-merge(nor_matrix,nor_data,by="gene")
}
rownames(nor_matrix)<-nor_matrix$gene
nor_matrix<-nor_matrix[,-1]

write.table(nor_matrix,"./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/si_ST_normal_matrix.txt",sep="\t")

normal_count<-nor_matrix
cell.data<-as.data.frame(colnames(normal_count))
cell.data$sample_info<-sapply(strsplit(cell.data[,1], "_"), function(x) paste(x[-length(x)], collapse = "_"))
colnames(cell.data)[1]<-"cell.names"


for (sam in c("RCC1","RCC2","RCC3","RCC4","RCC4f","RCC5")){
print(sam)

SCTdata<-read.table(paste0("./simulation/GSE224630_renal/si_ST/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
norm.cell.names<-cell.data$cell.names[which(cell.data$sample_info==sam)]
spatial_count<-SCTdata[,which(!(colnames(SCTdata) %in% norm.cell.names))]
spot_to_keep <- apply(spatial_count, 1, function(row) length(which(row > 0)) >= 100)
spatial_filtered <-spatial_count[spot_to_keep,]
rows_to_remove <- grepl("^MT-", rownames(spatial_filtered))
# Subset the dataframe to keep only the rows that don't start with "MT-"
spatial_filtered <- spatial_filtered[!rows_to_remove,]
print(dim(spatial_filtered))
spatial_filtered$gene<-rownames(spatial_filtered) 
normal_count$gene<-rownames(normal_count) 
all_count<-merge(spatial_filtered,normal_count,by="gene")
rownames(all_count)<-all_count$gene
all_count<-all_count[,-1]
print(length(intersect(colnames(all_count),cell.data$cell.names)))

print("start cnv")
cnv_data<-copykat(rawmat=all_count,id.type="S", LOW.DR=0.03, UP.DR=0.2,ngene.chr=5,sam.name=paste0(sam,"_st_0803"),KS.cut=0.1,win.size=25, norm.cell.names= cell.data$cell.names)
}

cnv_res<-data.frame()
for (sam in c("RCC1","RCC2","RCC3","RCC4","RCC4f","RCC5")){
print(sam)
cnv_data<-read.table(paste0("./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/cnv/",sam,"_st_0803_copykat_prediction.txt"),header=T,sep="\t")
cnv_data<-cnv_data%>%filter(copykat.pred!="not.defined"& grepl(paste0("^",sam,"_"), cell.names))
cnv_res<-rbind(cnv_res,cnv_data)
}
write.table(cnv_res,"./simulation/GSE224630_renal/si_ST/dif_Mspots_percet/cnv/5_sample_cnv_res.txt",sep="\t")

####################GSE181919_HNSC########################

scmeta<-read.table("./simulation/GSE181919_HNSC/GSE181919_Barcode_metadata.txt",header=T,sep="\t")
counts <- fread("./simulation/GSE181919_HNSC/GSE181919_UMI_counts.txt.gz")
counts2<-as.matrix(counts[,2:ncol(counts)])
rownames(counts2)<-as.character(counts$V1)
colnames(counts2)<-gsub("\\.", "-", colnames(counts2))
seurat_object <- CreateSeuratObject(counts = counts2)
seurat_object@meta.data$sample.id<-scmeta$sample.id
#seurat_object$orig.ident<-sapply(strsplit(rownames(seurat_object@meta.data), "-"),function(x) unlist(x)[2])
samples<-unique(scmeta$sample.id)[grepl("^C", unique(scmeta$sample.id))]

n_spots <- 1000
min_cells_per_spot <- 3
max_cells_per_spot <- 10
max_cell_types_per_spot <- 3


for (sam in  samples){
print(sam)
#M_spots_count <- round(M_appearance_rates[[sam]] * n_spots)
sc_object <-subset(seurat_object,subset =sample.id== sam)
#print(paste0(unique(seurat_object$orig.ident)[sam],",",unique(scmeta$sample)[sam]))
n_cells <- 10
test_scmeta<-scmeta%>%filter(sample.id==sam)
#test_scmeta<-test_scmeta%>%filter(Type!="Unresolved")
test_scmeta$Type<-test_scmeta$cell.type

rm_type<-as.data.frame(table(test_scmeta$Type))%>%filter(Freq<10)%>%.$Var1%>%as.character()
test_scmeta<-test_scmeta%>%filter(!Type %in% rm_type)

#test_scmeta<-test_scmeta[which(test_scmeta$nCount_RNA>1000),]
cell_types <-unique(test_scmeta$Type)
if(!"Malignant.cells" %in% cell_types){
next}

cell_data <-sc_object$RNA@data

# Create an empty list to store the simulated spots
spots <- vector("list", n_spots)
spots_celltype_info <- vector("list", n_spots)
# Simulate each spot
cell_type_table <- matrix(NA, nrow = n_spots, ncol = 6, 
                          dimnames = list(NULL, c("CellType1", "CellType2", "CellType3", "Prop1", "Prop2", "Prop3")))
set.seed(123)  
for (i in 1:n_spots) {
 # Randomly decide the number of cells in this spot
  n_cells_in_spot <- sample(min_cells_per_spot:max_cells_per_spot, 1)  
  # Randomly decide the number of cell types in this spot, not exceeding the max_cell_types_per_spot
  n_cell_types_in_spot <- sample(1:max_cell_types_per_spot, 1)
  #print(n_cell_types_in_spot)
  # Randomly decide the cell types in this spot
  cell_types_in_spot <- sample(cell_types, n_cell_types_in_spot, replace = FALSE)
  #print(paste0("cell_type:",length(unique(cell_types_in_spot))))
  # Generate cell type proportions using the Dirichlet distribution
  cell_type_proportions <- rdirichlet(1, rep(1, length(cell_types_in_spot)))[1,] 
  cell_type_table[i, 1:length(cell_types_in_spot)] <- cell_types_in_spot
  cell_type_table[i, (1:length(cell_type_proportions)) + 3] <- cell_type_proportions
  # For each cell type, sample the corresponding number of cells from the original data
  spot_cells <- lapply(1:n_cells_in_spot, function(x) {
    cell_type <- sample(cell_types_in_spot, 1, prob = cell_type_proportions)
    test_scmeta[sample(which(test_scmeta$Type == cell_type), 1), ]
  })
  spots_celltype_info[[i]] <- sapply(spot_cells, function(cell) cell$Type)
  # Combine the cells into a spot using the mean
  spots[[i]] <- apply(cell_data[,sapply(spot_cells, function(cell) rownames(cell))],1,mean)
}
rownames(cell_type_table)<-paste0("spot", 1:nrow(cell_type_table))
write.table(cell_type_table,paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/", sam,"si_ST_cell_type_propotion.txt"),sep="\t")
  
simulated_data <- do.call(rbind, spots)
rownames(simulated_data) <- paste0("spot", 1:nrow(simulated_data))
simulated_data<-t(simulated_data)
# Combine the cell type info into a data frame
spots_celltype_info <- lapply(spots_celltype_info, function(x) paste(x, collapse = ","))
# Convert the list into a dataframe
spots_celltype_info_df <- data.frame(Spot = 1:n_spots,CellTypes = unlist(spots_celltype_info))
spots_celltype_info_df$Num_Cell_Types <- sapply(spots_celltype_info, function(x) length(unique(strsplit(x, ",")[[1]])))
spots_celltype_info_df$Num_Cells <- sapply(spots_celltype_info, function(x) length(strsplit(x, ",")[[1]]))
rownames(spots_celltype_info_df) <- paste0("spot", 1:nrow(spots_celltype_info_df))
print(table(spots_celltype_info_df$Num_Cell_Types))
write.table(simulated_data,paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/", sam,"si_ST_mtx.txt"),sep="\t")
write.table(spots_celltype_info_df,paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/", sam,"si_ST_cell_type.txt"),sep="\t")
}

gsva_res<-data.frame()
for (sam in samples){
print(sam)
SCTdata<-read.table(paste0("./simulation/GSE181919_HNSC/si_ST/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
SCTdata<-as.matrix(SCTdata)
print(range(apply(SCTdata,2,function(x)length(which(x>0)))))
print(range(rowSums(SCTdata)))
print(range(colSums(SCTdata)))
gsva_score<-as.data.frame(t(gsva(SCTdata, markerlist, parallel.sz=1)))
gsva_res<-rbind(gsva_res,gsva_score)
write.table(gsva_res,"./simulation/GSE181919_HNSC/si_ST/gsva_htca_cancer_res_6.txt",sep="\t")
}
samples2<-samples[-which(samples=="C84")]
get_frequencies <- function(vector_string) {
  vector_elements <- strsplit(as.character(vector_string), split = ",")[[1]]
  return(table(vector_elements))
}

for (sam in samples2){
print(sam)
spotdata<-read.table(paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type.txt"),sep="\t",header=T)
SCTdata<-read.table(paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
frequencies <- lapply(spotdata$CellTypes, get_frequencies)
# Initialize a list to store the frequency data for each row
frequency_list <- list()
# Iterate through the frequencies and convert to data frame
for (i in seq_along(frequencies)) {
  frequency_list[[i]] <- as.data.frame(as.table(frequencies[[i]]))
}
frequency_df <- data.frame(element = character(), frequency = integer(),  percent=numeric(),row_name = character())
# Iterate through the frequency list, adding the required information to the final data frame
for (i in seq_along(frequency_list)){
  freq_table <- frequency_list[[i]]
  elements <- freq_table$vector_elements
  frequencies <- freq_table$Freq
  percentages <- frequencies / sum(frequencies) 
  row_df <- data.frame(
    element = elements,
    frequency = frequencies,
	percent=percentages,
    row_name = rep(rownames(spotdata)[i], length(elements)), # Repeat the row name for each element
    stringsAsFactors = FALSE
  )
  frequency_df <- rbind(frequency_df, row_df)
}
frequency_df$dataset<-sam
nor_spot<- frequency_df %>%
  group_by(row_name) %>%
  filter(!any(element == "Malignant.cells")) %>%
  select(row_name) %>%
  distinct()
set.seed(222)
nor_spots <- sample(nor_spot$row_name, 10, replace = FALSE)
nor_spots<-paste0(sam,"_",nor_spots)
nor_data<-SCTdata[,nor_spots]
write.table(frequency_df,paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/",sam,"si_ST_cell_type_summmary.txt"),sep="\t")
write.table(nor_data,paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/",sam,"si_ST_normal_matrix.txt"),sep="\t")
}
nor_matrix<-read.table(paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/",samples2[1],"si_ST_normal_matrix.txt"),sep="\t")
nor_matrix$gene<-rownames(nor_matrix)
for (sam in samples2[2:length(samples2)]){
print(sam)
nor_data<-read.table(paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/",sam,"si_ST_normal_matrix.txt"),sep="\t")
nor_data$gene<-rownames(nor_data)
nor_matrix<-merge(nor_matrix,nor_data,by="gene")
}
rownames(nor_matrix)<-nor_matrix$gene
nor_matrix<-nor_matrix[,-1]

write.table(nor_matrix,"./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/si_ST_normal_matrix.txt",sep="\t")

normal_count<-nor_matrix
cell.data<-as.data.frame(colnames(normal_count))
cell.data$sample_info<-sapply(strsplit(cell.data[,1], "_"), function(x) paste(x[-length(x)], collapse = "_"))
colnames(cell.data)[1]<-"cell.names"
for (sam in samples2[8:19]){
print(sam)
SCTdata<-read.table(paste0("./simulation/GSE181919_HNSC/si_ST/dif_Mspots_percet/",sam,"si_ST_mtx.txt"),sep="\t",header=T)
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
norm.cell.names<-cell.data$cell.names[which(cell.data$sample_info==sam)]
spatial_count<-SCTdata[,which(!(colnames(SCTdata) %in% norm.cell.names))]
gene_to_keep <- apply(spatial_count, 1, function(row) length(which(row > 0)) >= 100)
spatial_filtered <-spatial_count[gene_to_keep,]
rows_to_remove <- grepl("^MT-", rownames(spatial_filtered))
# Subset the dataframe to keep only the rows that don't start with "MT-"
spatial_filtered <- spatial_filtered[!rows_to_remove,]
print(dim(spatial_filtered))
spatial_filtered$gene<-rownames(spatial_filtered) 
normal_count$gene<-rownames(normal_count) 
all_count<-merge(spatial_filtered,normal_count,by="gene")
rownames(all_count)<-all_count$gene
all_count<-all_count[,-1]
print(length(intersect(colnames(all_count),cell.data$cell.names)))

print("start cnv")
cnv_data<-copykat(rawmat=all_count,id.type="S", LOW.DR=0.03, UP.DR=0.2,ngene.chr=5,sam.name=paste0(sam,"_st_1017_pan"),KS.cut=0.1,win.size=25, norm.cell.names= cell.data$cell.names)
}

cnv_res<-data.frame()
for (sam in samples){
print(sam)

cnv_data<-read.table(paste0("./simulation/GSE181919_HNSC/si_ST/cnv/",sam,"_st_0803_copykat_prediction.txt"),header=T,sep="\t")
cnv_data<-cnv_data%>%filter(copykat.pred!="not.defined"& grepl(paste0("^",sam), cell.names))
cnv_res<-rbind(cnv_res,cnv_data)
}
write.table(cnv_res,"./simulation/GSE181919_HNSC/si_ST/cnv/20_sample_cnv_res.txt",sep="\t")
