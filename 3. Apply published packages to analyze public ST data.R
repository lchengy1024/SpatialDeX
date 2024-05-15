########A breast cancer ST dataset (PMID:34493872) as an example#######

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(MCPcounter)
library(SPOTlight)
library(Matrix)
library(ggpubr)
library(png)
library(CARD)
library(mistyR)
library(stringr)
library(GSVA)
library(harmony)
library(stringr)
library(copykat)
library(reshape2)
library(SeuratDisk)
library(Giotto)
library(spacexr)
library(STdeconvolve)
library(SpaCET)

############tme marker and cancer marker for GSVA##################

markerlist<-read.table("./sc_pan_cancer/TME_HTCA_CANCER_marerklist.txt",sep="\t",header=T)
####remove markers appear in multiple cell types#####
markerlist2<-data.frame()
for (i in 1:ncol(markerlist)){
tmp<-as.data.frame(markerlist[,i])
markerlist2<-rbind(markerlist2,tmp)
}
a<-as.data.frame(table(markerlist2))
a<-arrange(a,Freq)
a<-a[which(a$Freq<3),]
markerlist<-as.list(markerlist)
markerlist<-lapply(markerlist,function(x) x<-x[which(x %in% a[,1])])
######gsva
file_folder<-list.files("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices")
file_name<-str_split_fixed(file_folder, fixed("_"),n=2)
file_name<-file_name[,1]
gsva_res<-data.frame()
for (i in 1:length(file_name)){
print(file_name[i])
st_object<-Read10X(paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",file_name[i],"_raw_feature_bc_matrix/"))
st_object <- CreateSeuratObject(st_object,min.cells = 3)
st_object[["percent.mt"]] <- PercentageFeatureSet(st_object, pattern = "^MT-")
st_object <- subset(st_object, subset = nFeature_RNA > 100  &  nCount_RNA>500& percent.mt < 20)
SCTdata<-as.matrix(st_object$RNA@data)
colnames(SCTdata)<-paste0(file_name[i],"_",colnames(SCTdata))
gsva_score<-as.data.frame(t(gsva(SCTdata, markerlist,verbose=FALSE, parallel.sz=1)))
gsva_res<-rbind(gsva_res,gsva_score)
write.table(gsva_res,"./spatial_public_data/34493872_BC/gsva_htca_cancer_res.txt",sep="\t")
}

####Giotto RCTD
test_sc<-Read10X("./spatial_public_data/34493872_BC/single_cell/sc_data")
scmeta<-read.csv("./spatial_public_data/34493872_BC/single_cell/Whole_miniatlas_meta.csv",header=T)
sam_name<-unique(scmeta$Patient)
file_name = c('CID4465','CID44971','CID4290A','CID4535')
for (i in 1:length(file_name)){
print(file_name[i])
print("running Giotto")
my_instructions <- createGiottoInstructions(save_plot=TRUE,
show_plot=TRUE,
return_plot=FALSE,
save_dir=paste0('./spatial_public_data/34493872_BC/output/giotto/minor/',file_name[i]))
st_object<-Read10X(paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",file_name[i],"_raw_feature_bc_matrix/"))
st_object <- CreateSeuratObject(st_object,min.cells = 3)
st_object[["percent.mt"]] <- PercentageFeatureSet(st_object, pattern = "^MT-")

st_object <- subset(st_object, subset = nFeature_RNA > 100  &  nCount_RNA>500& percent.mt < 20)
SCTdata<-as.matrix(st_object$RNA@data)
spatial_gitto<-createGiottoObject(raw_exprs = SCTdata,spatial_locs = NULL,
                                           instructions = my_instructions)									   
spatial_gitto <- filterGiotto(gobject = spatial_gitto,
                                     expression_threshold = 1,
                                     feat_det_in_min_cells = 10,
                                     min_det_feats_per_cell = 100,
                                     expression_values = c('raw'),
                                     verbose = T)
spatial_gitto <- normalizeGiotto(gobject = spatial_gitto)
spatial_gitto <- calculateHVF(gobject = spatial_gitto)	
gene_metadata <- fDataDT(spatial_gitto)
featgenes = gene_metadata[hvf == 'yes']$feat_ID
spatial_gitto <- runPCA(gobject = spatial_gitto, genes_to_use = featgenes, scale_unit = F)
spatial_gitto <- runUMAP(spatial_gitto, dimensions_to_use = 1:10)

spatial_gitto <- createNearestNetwork(gobject = spatial_gitto, dimensions_to_use = 1:10, k = 15)
spatial_gitto <- doLeidenCluster(gobject = spatial_gitto, resolution = 0.4, n_iterations = 1000)
print(colnames(pDataDT(spatial_gitto)))	
test_scmeta<-scmeta%>%filter(Patient==file_name[i])
test_scmeta<-test_scmeta[-1,]
sc_meta<-test_scmeta[which(test_scmeta$nCount_RNA>1000),]
rownames(sc_meta)<-sc_meta$NAME
sc_count<-test_sc[,sc_meta$NAME]
# gene_rm<-data.frame(table(rownames(sc_count)))%>%filter(Freq>=2)%>%.$Var1%>%as.character()
# sc_count_ref<-sc_count[-which(rownames(sc_count)%in% gene_rm),sc_meta$cell_id]
sc_gitto <- createGiottoObject(raw_exprs = sc_count,instructions = my_instructions)
sc_gitto <- addCellMetadata(sc_gitto,new_metadata = sc_meta)
cell_count<-as.data.frame(table(pDataDT(sc_gitto)$celltype_minor))%>%filter(Freq<10)%>%.$Var1%>%as.character()
cell_id <- pDataDT(sc_gitto)%>%filter(!(celltype_minor %in% cell_count))%>%.$cell_ID
sc_gitto <- subsetGiotto(sc_gitto, cell_ids = cell_id )
sc_gitto <- normalizeGiotto(gobject = sc_gitto)
markers_scran = findMarkers_one_vs_all(gobject=sc_gitto, method="scran",
                                    expression_values="normalized", cluster_column='celltype_minor', min_feats=3)
top_markers <- markers_scran[, head(.SD, 15), by="cluster"]
print(length(top_markers))
DWLS_matrix<-makeSignMatrixDWLSfromMatrix(matrix = as.matrix(slot(sc_gitto@expression$cell$rna$normalized, 'exprMat')),
                                        cell_type = pDataDT(sc_gitto)$celltype_minor,
                                        sign_gene = top_markers$feats)
print(colnames(pDataDT(spatial_gitto)))
spatial_gitto = runDWLSDeconv(gobject = spatial_gitto, expression_values = c("normalized"),sign_matrix = DWLS_matrix)
gitto_res <- as.data.frame(spatial_gitto@spatial_enrichment$cell$rna$DWLS@enrichDT)
gitto_res$cell.names<-paste0(file_name[i],"_",rownames(gitto_res))
write.table(gitto_res,paste0("./spatial_public_data/34493872_BC/output/giotto/minor/",file_name[i],"_giotto_res.txt"),sep="\t")
##########RCTD
print("running rctd")
st_counts<-SCTdata
coords<-as.data.frame(matrix(1,nrow=dim(SCTdata)[2],ncol=2))
colnames(coords)<-c("x","y");rownames(coords)<-colnames(SCTdata)
nUMI <- st_object$nCount_RNA; names(nUMI) <- rownames(st_object@meta.data)
st_data <- SpatialRNA(coords, st_counts, nUMI)
#barcodes <- colnames(st_data@counts)
cell_types <- sc_meta$celltype_minor;names(cell_types) <-sc_meta$NAME # convert to factor data type
nUMI_sc <- as.numeric(sc_meta$nCount_RNA); names(nUMI_sc) <- sc_meta$NAME # create nUMI named list
### Create the Reference object
sc_count2<-sc_count[,cell_id]
nUMI_sc<-nUMI_sc[which(names(nUMI_sc)%in%cell_id)]
cell_types<-as.factor(cell_types[which(names(cell_types)%in%cell_id)])
reference <- Reference(sc_count2, cell_types, nUMI_sc)
myRCTD <- create.RCTD(st_data , reference, max_cores = 1,CELL_MIN_INSTANCE = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
results <- myRCTD@results
norm_weights = normalize_weights(results$weights)
norm_weights<-norm_weights[,sort(colnames(norm_weights))]
write.table(norm_weights,paste0("./spatial_public_data/34493872_BC/output/RCTD/minor/",file_name[i],"_full_mode_res_minor.txt"),sep="\t")
rank_res<-data.frame()
for (n in 1:nrow(norm_weights)){
tmp<-norm_weights[n,]
tmp<-sort(tmp,decreasing = TRUE)
rank_res[n,1:3]<-tmp[1:3]
rank_res[n,4:6]<-names(tmp)[1:3]
}
rownames(rank_res)<-rownames(norm_weights)
colnames(rank_res)<-c("sc_score1","sc_score2","sc_score3","sc_type1","sc_type2","sc_type3")
rank_res$cell.names<-paste0(file_name[i],"_",rownames(rank_res))
write.table(rank_res,paste0("./spatial_public_data/34493872_BC/output/RCTD/minor/",file_name[i],"_rctd_rank_res_0607.txt"),sep="\t")
}

########CARD##############
########reference-based mode#########
test_sc<-Read10X("./spatial_public_data/34493872_BC/single_cell/sc_data")
scmeta<-read.csv("./spatial_public_data/34493872_BC/single_cell/Whole_miniatlas_meta.csv",header=T)
file_name = c('CID4465','CID44971','CID4290','CID4535')

for (sam in file_name){
print(sam)
st_object<-Read10X(paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",sam,"_raw_feature_bc_matrix/"))
st_object <- CreateSeuratObject(st_object,min.cells = 3)
st_image<-Read10X_Image(
  paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",sam,"_raw_feature_bc_matrix/spatial")
)
Seurat::DefaultAssay(object = st_image) <- 'RNA'
st_object<- subset(st_object, cells = Cells(st_image))
st_object[["image"]] <- st_image
st_object[["percent.mt"]] <- PercentageFeatureSet(st_object, pattern = "^MT-")
# st_object <- subset(st_object, subset = nFeature_RNA > 100  &  nCount_RNA>500& percent.mt < 20)
st_object@meta.data$cell.names<-paste0(sam,"_",rownames(st_object@meta.data))
SCTdata<-as.matrix(st_object$RNA@data)
spatial_location<-data.frame(st_object@images$image@coordinates)
spatial_location<-spatial_location%>%filter(tissue==1)
spatial_location<-spatial_location[,c(2,3)]
colnames(spatial_location)<-c("x","y")
if(sam=="CID4290"){
test_scmeta<-scmeta%>%filter(Patient==paste0(sam,"A"))
}else{
test_scmeta<-scmeta%>%filter(Patient==sam)
}
test_scmeta<-test_scmeta[-1,]
sc_meta<-test_scmeta[which(test_scmeta$nCount_RNA>1000),]
rownames(sc_meta)<-sc_meta$NAME
sc_count<-test_sc[,sc_meta$NAME]
meta_data<-sc_meta
meta_data<-meta_data%>%filter(celltype_major!="Unresolved")
rm_type<-as.data.frame(table(meta_data$celltype_major))%>%filter(Freq<10)%>%.$Var1%>%as.character()
meta_data<-meta_data%>%filter(!celltype_major %in% rm_type)
sc_data <-sc_count[,meta_data$NAME]

CARD_obj = createCARDObject(
	sc_count = sc_data,
	sc_meta = meta_data,
	spatial_count = SCTdata,
	spatial_location = spatial_location,
	ct.varname = "celltype_major",
	ct.select = unique(meta_data$celltype_major),
	sample.varname = "Patient",
	minCountGene = 100,
	minCountSpot = 5) 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)	
res<-CARD_obj@Proportion_CARD
write.table(res,paste0("./spatial_public_data/34493872_BC/output/CARD/reference_based_major/",sam,"_CARD_res.txt"),sep="\t")
}
########reference free mode##########
for (sam in sample_in){
st_object<-Read10X(paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",sam,"_raw_feature_bc_matrix/"))
st_object <- CreateSeuratObject(st_object,min.cells = 3)
st_image<-Read10X_Image(
  paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",sam,"_raw_feature_bc_matrix/spatial")
)
Seurat::DefaultAssay(object = st_image) <- 'RNA'
st_object<- subset(st_object, cells = Cells(st_image))
st_object[["image"]] <- st_image
st_object <- subset(st_object, subset = nFeature_RNA > 100  &  nCount_RNA>500)
spatial_coord <- data.frame(st_object@images$image@coordinates)%>%
        tibble::rownames_to_column("barcodeID") %>% dplyr::mutate(imagerow_scaled = imagerow *
       st_object@images$image@scale.factors$lowres, imagecol_scaled = imagecol *
        st_object@images$image@scale.factors$lowres)
coord_data<-spatial_coord[,c(1,7,8)]
rownames(coord_data)<-coord_data[,1]
coord_data<-coord_data[,-1]
colnames(coord_data)<-c("x","y")

spatial_count<-st_object@assays$RNA@counts 
spatial_location<-data.frame(coord_data)

CARDfree_obj = createCARDfreeObject(
	markerList = markerlist,
	spatial_count = spatial_count,
	spatial_location = spatial_location,
	minCountGene = 100,
	minCountSpot =5) 
CARDfree_obj = CARD_refFree(CARDfree_obj)	
res<-CARDfree_obj@Proportion_CARD	
#saveRDS(CARDfree_obj, file = paste0("./spatial_public_data/test/CARD_obj_merge_T_tme_",cancer_types[i],".rds"))
write.table(res,paste0("./spatial_public_data/34493872_BC/CARD/",sam,"_CARD_res_1129.txt"),sep="\t")
}
                   
##########SpaCET############
data.merge<-readRDS("./spatial_public_data/34493872_BC/spatial_Seurat_merge_allBC2.rds")
sample_in<-c("CID4535","CID4465","CID44971","CID4290","1142243F","1160920F")
for (sam in sample_in){
st_object<-Read10X(paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",sam,"_raw_feature_bc_matrix/"))
st_object <- CreateSeuratObject(st_object,min.cells = 3)
st_image<-Read10X_Image(
  paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",sam,"_raw_feature_bc_matrix/spatial")
)
Seurat::DefaultAssay(object = st_image) <- 'RNA'
st_object<- subset(st_object, cells = Cells(st_image))
st_object[["image"]] <- st_image
st_object <- subset(st_object, subset = nFeature_RNA > 100  &  nCount_RNA>500)
spatial_coord <- data.frame(st_object@images$image@coordinates)%>%
        tibble::rownames_to_column("barcodeID") %>% dplyr::mutate(imagerow_scaled = imagerow *
       st_object@images$image@scale.factors$lowres, imagecol_scaled = imagecol *
        st_object@images$image@scale.factors$lowres)
coord_data<-spatial_coord[,c(1,7,8)]
rownames(coord_data)<-coord_data[,1]
coord_data<-coord_data[,-1]

spatial_count<-st_object@assays$RNA@counts 
SpaCET_obj <- create.SpaCET.object(
  counts=spatial_count,
  spotCoordinates=coord_data,
  imagePath=NA,
   platform = "Visium"
)
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType="PANCAN", coreNo=8)
res<-data.frame(t(SpaCET_obj@results$deconvolution$propMat))
write.table(res,paste0("./comparison/SpaCET/res_",sam,".txt"),sep="\t")
}

##########STdeconvolve#############
test_sc<-Read10X("./spatial_public_data/34493872_BC/single_cell/sc_data")
scmeta<-read.csv("./spatial_public_data/34493872_BC/single_cell/Whole_miniatlas_meta.csv",header=T)
file_name = c('CID4465','CID44971','CID4290','CID4535')
file_name = c('CID4290','CID4535')

for (sam in file_name){
print(sam)
st_object<-Read10X(paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",sam,"_raw_feature_bc_matrix/"))
st_object <- CreateSeuratObject(st_object,min.cells = 3)
st_image<-Read10X_Image(
  paste0("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices/",sam,"_raw_feature_bc_matrix/spatial")
)
Seurat::DefaultAssay(object = st_image) <- 'RNA'
st_object<- subset(st_object, cells = Cells(st_image))
st_object[["image"]] <- st_image

st_object[["percent.mt"]] <- PercentageFeatureSet(st_object, pattern = "^MT-")
# st_object <- subset(st_object, subset = nFeature_RNA > 100  &  nCount_RNA>500& percent.mt < 20)
st_object@meta.data$cell.names<-paste0(sam,"_",rownames(st_object@meta.data))
SCTdata<-as.matrix(st_object$RNA@data)

if(sam=="CID4290"){
test_scmeta<-scmeta%>%filter(Patient==paste0(sam,"A"))
}else{
test_scmeta<-scmeta%>%filter(Patient==sam)
}
test_scmeta<-test_scmeta[-1,]
sc_meta<-test_scmeta[which(test_scmeta$nCount_RNA>1000),]
rownames(sc_meta)<-sc_meta$NAME
sc_count<-test_sc[,sc_meta$NAME]
meta_data<-sc_meta
meta_data<-meta_data%>%filter(celltype_major!="Unresolved")
rm_type<-as.data.frame(table(meta_data$celltype_major))%>%filter(Freq<10)%>%.$Var1%>%as.character()
meta_data<-meta_data%>%filter(!celltype_major %in% rm_type)
sc_data <-sc_count[,meta_data$NAME]
aver_all<-data.frame(rownames(sc_data))
for (ctype in unique(meta_data$celltype_major)){
cell_b<-meta_data%>%filter(celltype_major==ctype)%>%rownames()
aver<-data.frame(apply(sc_data[,cell_b],1,mean))
colnames(aver)<-ctype
aver_all<-cbind(aver_all,aver)
}
#rm_gene<-as.data.frame(table(aver_all[,1]))%>%filter(Freq>1)%>%.$Var1
#aver_all<-aver_all[-which(aver_all[,1]%in% rm_gene),]
rownames(aver_all)<-aver_all[,1]
aver_all<-aver_all[,-1]
colnames(SCTdata)<-paste0(sam,"_",colnames(SCTdata))
SCTdata<-as.matrix(SCTdata)
SCTdata<- round(SCTdata, digits = 0)
corpus <- restrictCorpus(SCTdata, removeAbove=1.0, removeBelow = 0.1)
all_zero_rows <- apply(corpus, 2, function(x) all(x == 0))
corpus_filtered <- corpus[,!all_zero_rows]
ldas <- fitLDA(t(as.matrix(corpus_filtered )), Ks = seq(2, 6, by = 1))
optLDA <- optimalModel(models = ldas, opt = "min")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
aver_all2<-aver_all[colnames(deconGexp),]
corMtx_beta <- getCorrMtx(m1 = as.matrix(deconGexp), # the deconvolved cell-type `beta` (celltypes x genes)
                          m2 = t(as.matrix(aver_all2)), # the reference `beta` (celltypes x genes)
                          type = "b") # "b" = comparing beta matrices, "t" for thetas
corMtx_beta<-data.frame(corMtx_beta)
corMtx_beta$max<-apply(corMtx_beta,1,function(x) names(x)[which.max(x)])
colnames(deconProp)<-corMtx_beta$max
write.table(deconProp,paste0("./spatial_public_data/34493872_BC/output/STdeconvolove/major/",sam,"_STdecon_res.txt"),sep="\t")
}

file_name = c('CID4465','CID44971','CID4290','CID4535')
pa<-read.table("./spatial_public_data/34493872_BC/output/34493872_BC_pathology.txt",sep="\t",header=T)
for (sam in file_name){
print(sam)
STdecon_res<-read.table(paste0("./spatial_public_data/34493872_BC/output/STdeconvolove/major/",sam,"_STdecon_res.txt"),check.names=FALSE,sep="\t")
ct_freq<-data.frame(table(colnames(STdecon_res)))
n_col<-ncol(STdecon_res)
rm_col<-which(colnames(STdecon_res)[1:n_col]%in%ct_freq$Var1[which(ct_freq$Freq>1)])

sum_list<-list()
li=1
for(nr in 1:nrow(ct_freq)){
if(ct_freq$Freq[nr]>1){
sum_list[[li]]<-apply(STdecon_res,1,function(x)sum(x[which(names(x)==ct_freq$Var1[nr])]))
names(sum_list)[li]<-as.character(ct_freq$Var1[nr])
li=li+1
}}
STdecon_res<-cbind(STdecon_res,data.frame(sum_list))
STdecon_res<-STdecon_res[,-rm_col]
STdecon_res$max<-apply(STdecon_res,1,function(x) names(x)[which.max(x)])
#write.table(STdecon_res,paste0("./spatial_public_data/34493872_BC/output/STdeconvolove/major/",sam,"_STdecon_max_res.txt"),sep="\t")
pa_tmp<-pa%>%filter(patientid==sam)
pa_tmp$ID<-paste0(sam,"_",pa_tmp$ID)
STdecon_res$ID<-rownames(STdecon_res)
STdecon_res2<-merge(STdecon_res,pa_tmp,by="ID")
STdecon_res2$max <- 
               ifelse(STdecon_res2$max == "Plasmablasts", "lymphocyte",
			   ifelse(STdecon_res2$max == "CAFs", "stroma",
                      ifelse(STdecon_res2$max == "Myeloid", "other", STdecon_res2$max)))					  
write.table(STdecon_res2,paste0("./spatial_public_data/34493872_BC/output/STdeconvolove/major/",sam,"_STdecon_max_res.txt"),sep="\t")
}
