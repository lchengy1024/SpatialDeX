###########functions###################
.libPaths(new = c("/home/xinyiliu/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))

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
library(infercnv)
library(mistyR)

library(GSVA)
library(harmony)
library(stringr)
library(copykat)

library(progeny)
library(philentropy)
library(reshape2)
library(SeuratDisk)
library(Giotto)
library(spacexr)
library(data.table)
library(Matrix)
library(xgboost)
library(caret)

library(plyr)
library(tidymodels)
library(themis)
library(vip)
library(stringr)
library(jaccard)
library(pROC)
library(mltools)
getsolvedmatrix <- function(Xm, Ym, lambda = 0.01) {
    # Amat=solve(Xmat,Ymat) # solve AX=B, or Xmat * A =Ymat
    TMmat_g = (t(Xm) %*% Xm) + lambda * diag(ncol(Xm))
    
    Amat_g = solve(TMmat_g) %*% t(Xm) %*% Ym
    return(Amat_g)
}
getsolvedmatrix_with_permutation_cell_label <- function(Xm, Ym, lambda = 0.01, npermutation = 10000) {
    set.seed(000)
	Amat_ret = getsolvedmatrix(Xm, Ym, lambda = lambda)
    Amat_ret_higher = matrix(rep(0, ncol(Amat_ret) * nrow(Amat_ret)), nrow = nrow(Amat_ret))
    rownames(Amat_ret_higher) = rownames(Amat_ret)
    colnames(Amat_ret_higher) = colnames(Amat_ret)
    # permute N times randomly shuffle cell labels
    for (npm in 1:npermutation) {
        if (npm%%5000 == 0) {
            message(paste("Permutation:", npm, "/", npermutation, "..."))
        }
        cells_shu = sample(rownames(Ym), nrow(Ym))
        Xm_s = Xm[cells_shu, ]
        Ym_s = Ym  # [cells_shu,]
        rownames(Ym_s) = cells_shu
        Amat_random = getsolvedmatrix(Xm_s, Ym_s, lambda = lambda)
        
        Amat_ret_higher = Amat_ret_higher + (abs(Amat_random) > abs(Amat_ret)) * 1
        # browser()
    }
    Amat_ret_higher = Amat_ret_higher/npermutation
    return(list(Amat_ret, Amat_ret_higher))
}

get_top_n_and_names <- function(row, n) {
  top_indices <- order(row, decreasing = TRUE)[1:n]
  bottom_indices <- order(row)[1:n]
  top_colnames <- colnames(Amat_s2)[top_indices]
  bottom_colnames <- colnames(Amat_s2)[bottom_indices]
  #c(top_colnames,bottom_colnames)
  c(top_colnames)
}
update_B_based_on_A <- function(B, A) {
  B_new <-  matrix(rep(0, nrow(B) * ncol(B)), nrow = nrow(B))
  rownames(B_new) = rownames(B)
  colnames(B_new) = colnames(B) 
  for (i in 1:nrow(A)) {
    row_name <- rownames(A)[i]
    for (j in 1:ncol(A)) {
      col_name <- as.character(A[i, j])
      B_new[row_name, col_name] <-  B[row_name, col_name]
    }
  } 
  return(B_new)
}

scale_data<-function(dat){
dat_new<-data.frame()
for (sam in unique(dat$sample_info)){
tmp<-dat%>%filter(sample_info==sam)
tmp[is.na(tmp)]<-0
tmp[,-ncol(dat)]<-apply(tmp[,-ncol(dat)],c(1,2),as.numeric)
tmp[,-ncol(dat)]<-apply(tmp[,-ncol(dat)],2,function(x) (x-min(x))/(max(x)-min(x)))
tmp[is.na(tmp)]<-0
dat_new<-rbind(dat_new,tmp)}
return(dat_new)}


createStratifiedFolds <- function(samples, n = 4) {
  # Define number of samples to be selected from each part for a fold
  sizes <- c(2, 2, 1, 1)
  folds <- vector("list", n)
  
  # Shuffle each part once
  shuffled1 <- sample(samples[1:8])
  shuffled2 <- sample(samples[9:16])
  shuffled3 <- sample(samples[17:20])
  shuffled4 <- sample(samples[21:24])
  
  # Split each shuffled part into 'n' folds
  split1 <- split(shuffled1, cut(seq_along(shuffled1), n, labels = FALSE))
  split2 <- split(shuffled2, cut(seq_along(shuffled2), n, labels = FALSE))
  split3 <- split(shuffled3, cut(seq_along(shuffled3), n, labels = FALSE))
  split4 <- split(shuffled4, cut(seq_along(shuffled4), n, labels = FALSE))
  
  # Combine the samples from each split to create folds
  for (i in 1:n) {
    folds[[i]] <- c(split1[[i]], split2[[i]], split3[[i]], split4[[i]])
  }
  
  return(folds)
}
                        
#################simulated training and validation data###############################
                        
cnv_res<-read.table("./simulation/pan_cancer/cnv/pan_sample_cnv_res.txt",sep="\t",header=T)
gsva1<-read.table("./simulation/GSE182227_OPSCC/si_ST/gsva_htca_cancer_res.txt",sep="\t")
gsva2<-read.table("./simulation/GSE165897_OVCA/si_ST/gsva_htca_cancer_res.txt",sep="\t")
gsva3<-read.table("./simulation/GSE181919_HNSC/si_ST/gsva_htca_cancer_res.txt",sep="\t")
gsva4<-read.table("./simulation/GSE224630_renal/si_ST/gsva_htca_cancer_res.txt",sep="\t")
gsva_res<-rbind(gsva1,gsva2)
gsva_res<-rbind(gsva_res,gsva3)
gsva_res<-rbind(gsva_res,gsva4)
gsva_res$cell.names<-rownames(gsva_res)
com_res<-merge(cnv_res,gsva_res,by="cell.names")
com_res$sample_info<-sapply(strsplit(com_res$cell.names, "_"), function(x) paste(x[-length(x)], collapse = "_"))
#all_com2<-read.table("./simulation/pan_cancer/pan_sample_cell_type_res.txt",sep="\t",header=T) 
all_com2<-read.table("./simulation/pan_cancer/pan_sample_cell_type_summary.txt",sep="\t",header=T) 

cell_map<-read.table("./simulation/pan_cancer/cell_type_map_0803.txt",sep="\t",header=T)
head(cell_map)
all_com2$CellType1 <- all_com2 %>%
  left_join(cell_map, by = c("CellType1" = "cell_label")) %>%
  .$uni_label

all_com2$CellType2 <- all_com2 %>%
  left_join(cell_map, by = c("CellType2" = "cell_label")) %>%
  .$uni_label

all_com2$CellType3 <- all_com2 %>%
  left_join(cell_map, by = c("CellType3" = "cell_label")) %>%
  .$uni_label
cell_type <-unique( unlist(lapply(all_com2[,c(1:3)], as.character)))
cell_type <- cell_type[!is.na(cell_type)]

ind_matrix = matrix(rep(FALSE, length(rownames(all_com2)) * length(cell_type)), nrow = length(rownames(all_com2)))
rownames(ind_matrix) = rownames(all_com2)
colnames(ind_matrix) = cell_type
for (i in colnames(ind_matrix)){
print(i)
for (j in 1:nrow(ind_matrix)){
num_i<-length(which(all_com2[j,]==i))
tmp_s<-as.numeric(all_com2[j,which(all_com2[j,]==i)+3])
if(i=="Malignant"){
if(num_i!=0 & length(which(tmp_s>0.1))==1){
ind_matrix[j,i]<-1
}else{
ind_matrix[j,i]<-0
}}
else{if(num_i!=0 & length(which(tmp_s>0.1))==1){
ind_matrix[j,i]<-1
}else{
ind_matrix[j,i]<-0
}
}
}}

a<-rowSums(ind_matrix==1)
sp<-names(a)[which(a!=0)]
ind_matrix2<-data.frame(ind_matrix[sp,])
ind_matrix2 <- ind_matrix2[!(ind_matrix2[,9] == 1 |ind_matrix2[,11] == 1 |ind_matrix2[,12] == 1 |ind_matrix2[,13] == 1 | ind_matrix2[,14] == 1 | ind_matrix2[,15] == 1), ]
b<-colSums(ind_matrix2==1)
ct<-names(b)[which(b!=0)]
ind_matrix2<-ind_matrix2[,which(colnames(ind_matrix2)%in%ct)]
ind_matrix2$sum<-rowSums(ind_matrix2==1)
ind_matrix2$class<-ifelse(ind_matrix2$Malignant==1 &ind_matrix2$sum==1,"tumor","other")
ind_matrix2$class<-ifelse(ind_matrix2$Malignant==1 &ind_matrix2$sum>1,"tumor_tme",ind_matrix2$class)
ind_matrix2$class<-ifelse(ind_matrix2$Malignant==0, "tme",ind_matrix2$class)
ind_matrix2$cell.names<-rownames(ind_matrix2)
ind_matrix2$sample_info<-sapply(strsplit(rownames(ind_matrix2), "_"), function(x) paste(x[-length(x)], collapse = "_"))

all_data<-merge(com_res,ind_matrix2,by="cell.names")
table(all_data$class,all_data$copykat.pred)
sample_all<-unique(all_data$sample_info.x)


n_permutation = 10000
set.seed(121212)
#sct<- c(0,1)
sct<- c(0)
combinations <- expand.grid(replicate(9, sct, simplify = FALSE))
LAMBDA=0.01
p=0.0001
lambda=1
data_oc<-all_data[,c(1,70:78)]
data_fea<-all_data[,c(1:54,69)]
data_fea$copykat.pred2<-ifelse(data_fea$copykat.pred=="aneuploid",1,0)
data_fea[is.na(data_fea)]<-0


data_fea$tme_av0<-apply(data_fea[,11:54],1,function(x) mean(x[which(x>0)]))
data_fea$tme_av0[is.na(data_fea$tme_av0)]<-0
data_fea$tme_count<-apply(data_fea[,11:54],1,function(x) length(which(x>0)))
data_fea$tme_av<-apply(data_fea[,11:54],1,mean)
data_fea$htca_av0<-apply(data_fea[,4:9],1,function(x) mean(x[which(x>0)]))
data_fea$htca_av0[is.na(data_fea$htca_av0)]<-0
data_fea$htca_count<-apply(data_fea[,4:9],1,function(x) length(which(x>0)))
data_fea$htca_av<-apply(data_fea[,4:9],1,mean)

data_fea2<-data.frame()
for (sam in unique(data_fea$sample_info.x)){
tmp<-data_fea%>%filter(sample_info.x==sam)
tmp[is.na(tmp)]<-0
tmp[,c(3:54,57,59,60,62)]<-apply(tmp[,c(3:54,57,59,60,62)],c(1,2),as.numeric)
tmp[,c(3:54,57,59,60,62)]<-apply(tmp[,c(3:54,57,59,60,62)],2,function(x) (x-min(x))/(max(x)-min(x)))
tmp[is.na(tmp)]<-0
data_fea2<-rbind(data_fea2,tmp)}
data_fea2$class<-all_data$class

rownames(data_oc)<-data_oc$cell.names
data_oc2<-data_oc[,-1]
# rms<-which(rowSums(data_oc2==1)==0|rowSums(data_oc2==1)>3)
# data_oc2<-data_oc2[-rms,]
data_oc2$sample_info<-sapply(strsplit(rownames(data_oc2), "_"), function(x) paste(x[-length(x)], collapse = "_"))


# data_c<-data_fea2
# rownames(data_c)<-data_c$cell.names
# data_c<-data_c[,-1]
# data_c<-data_c[rownames(data_oc2),]
# data_c$class5<-ifelse(data_c$class=="tumor","tumor","non_tumor")

multi_res<-data.frame()
mm=1
multi_res_a<-data.frame()
mma=1
res_sample<-data.frame()
res<-data.frame()
mlp2_res<-data.frame()
ll2=1
auc_res<-data.frame()
n=1
mean_auc<-data.frame()
ll=1
j_res<-data.frame(c(1:9))
percent<- c(0.1,0.2)
# generate all possible combinations of values
set.seed(131313)
percent_com <- expand.grid(replicate(9, percent, simplify = FALSE))
sam_acc_res<-data.frame()

cv_i=4
auc_res_sam<-data.frame()
print(paste0("sample combanition: ",cv_i))
sample_all<-unique(all_data$sample_info.x)
set.seed(321321321)
part1 <- sample_all[1:19]
part2 <- sample_all[20:38]
part3 <- sample_all[39:43]
part4 <- sample_all[44:49]

# Initialize a list to store results
results <- list()

# Repeat the sampling process 3 times
for (i in 1:10) {
  # Randomly select names from each part
  selected1 <- sample(part1, 8)
  selected2 <- sample(part2, 8)
  selected3 <- sample(part3, 4)
  selected4 <- sample(part4, 4) # From the last 10
  
  # Combine selected names
  combined <- c(selected1, selected2, selected3, selected4)
  
  # Store in the results list
  results[[i]] <- combined
}

# View the results
results
set.seed(3322211)
cv_results <- createStratifiedFolds(results[[cv_i]])
sam_t<-unlist(cv_results)
sam_v<-sample_all[-which(sample_all %in% sam_t)]
#sam_v<-sam_v[-which(sam_v %in% sample_to_remove)]
print(sam_t)
print(sam_v)

Y_mat_all<-data_fea2[,c(1,3:54,56,55)]
rownames(Y_mat_all)<-Y_mat_all$cell.names
Y_mat_all<-Y_mat_all[,-1]
Y_mat_all2<-Y_mat_all[rownames(data_oc2),]

ind_matrix_t1<-data_oc2%>%filter(sample_info %in% sam_t)
Y_mat_all_t1<-Y_mat_all2%>%filter(sample_info.x %in% sam_t)


Y_mat_all_v1<-Y_mat_all2%>%filter(sample_info.x %in% sam_v)
ind_matrix_v1<-data_oc2%>%filter(sample_info %in% sam_v)

ind_matrix_t1<-ind_matrix_t1[,-10]
Y_mat_all_t1<-Y_mat_all_t1[,-54]

print("####################start matrix decomposition###############")

ind_matrix_t<-as.matrix(ind_matrix_t1)
print(colSums(ind_matrix_t==1))
Y_mat_all_t<-as.matrix(Y_mat_all_t1)
Y_mat_all_t<-apply(Y_mat_all_t,c(1,2),as.numeric)
Yt<-Y_mat_all_t
set.seed(112233)
Amat_pm_lst = getsolvedmatrix_with_permutation_cell_label(ind_matrix_t, Y_mat_all_t, lambda = LAMBDA, npermutation = n_permutation)
Amat_s = Amat_pm_lst[[1]]
Amat_pval = Amat_pm_lst[[2]]
Amat_s<-t(scale(t(Amat_s),center=TRUE,scale=TRUE))
Amat_s[is.na(Amat_s)]<-0
rcom=1
Amat_s2 = matrix(rep(0, nrow(Amat_s) * ncol(Amat_s)), nrow = nrow(Amat_s))
rownames(Amat_s2) = rownames(Amat_s)
colnames(Amat_s2) = colnames(Amat_s) 
for (rn in 1:nrow(Amat_s)){
for (cn in 1:ncol(Amat_s)){
if(Amat_pval[rn,cn]<p & abs(Amat_s[rn,cn])>combinations[rcom,rn]){
Amat_s2[rn,cn]<-Amat_s[rn,cn]}
else{Amat_s2[rn,cn]<-Amat_s2[rn,cn]}
}
}
A<-Amat_s2
set.seed(123123)
Xt<-Yt%*%t(A)%*%solve((A %*% t(A)) +lambda * diag(nrow(A)))
Xt<-as.data.frame(Xt)
Xt[Xt<0]=0
Xt2<-Xt


ind_matrix_v1<-ind_matrix_v1[,-10]
ind_matrix_v<-as.matrix(ind_matrix_v1)
print(colSums(ind_matrix_v==1))

Y_mat_all_v1<-Y_mat_all_v1[,-54]
Y_mat_all_v<-as.matrix(Y_mat_all_v1)

Y_mat_all_v<-apply(Y_mat_all_v,c(1,2),as.numeric)
Y<-Y_mat_all_v
X<-Y%*%t(A)%*%solve((A %*% t(A)) +lambda * diag(nrow(A)))
X<-as.data.frame(X)
X[X<0]=0

X2=X
replace_smallest_with_zero <- function(row) {
  smallest_indices <- order(row)[1:6]  # Indices of the 5 smallest values
  #smallest_indices <- which(row<0.1)
  row[smallest_indices] <- 0          # Replace with 0
  return(row)
}
X2 <- as.data.frame(t(apply(X2, 1, replace_smallest_with_zero)))
Xt2 <- as.data.frame(t(apply(Xt2, 1, replace_smallest_with_zero)))
X2<-t(apply(X2,1,function(x) x/sum(x)))
Xt2<-t(apply(Xt2,1,function(x) x/sum(x)))
X2[is.nan(X2)] <- 0
Xt2[is.nan(Xt2)] <- 0
#################################

#combined_cell_propo<-read.table("./simulation/pan_cancer/pan_sample_cell_type_res.txt",sep="\t",header=T)  
combined_cell_propo<-read.table("./simulation/pan_cancer/pan_sample_cell_type_summary.txt",sep="\t")
          
cell_map<-read.table("./simulation/pan_cancer/cell_type_map_0803.txt",sep="\t",header=T)
combined_cell_propo[is.na(combined_cell_propo)] <- 0

combined_cell_propo$CellType1 <- combined_cell_propo %>%
  left_join(cell_map, by = c("CellType1" = "cell_label")) %>%
  .$uni_label
combined_cell_propo$CellType2 <- combined_cell_propo %>%
  left_join(cell_map, by = c("CellType2" = "cell_label")) %>%
  .$uni_label

combined_cell_propo$CellType3 <-combined_cell_propo %>%
  left_join(cell_map, by = c("CellType3" = "cell_label")) %>%
  .$uni_label
combined_cell_propo<-as.matrix(combined_cell_propo)

top1 <- apply(combined_cell_propo, 1, function(x) {
  name_cols <- x[1:3]
  max_val <- which.max(x[4:6])
  return (name_cols[max_val])})

top2 <- apply(combined_cell_propo, 1, function(x) {
  names_1_to_3 <- x[1:3]
  values_4_to_6 <- as.numeric(x)[4:6]
  
sorted_values <- sort(values_4_to_6, decreasing = TRUE)
second_largest <- sorted_values[2]
if (second_largest< 0.3) {
    return(NA)}else{

index<-which(values_4_to_6==second_largest)[1]
return(names_1_to_3[index])}

})

top3 <- apply(combined_cell_propo, 1, function(x) {
  names_1_to_3 <- x[1:3]
  values_4_to_6 <- as.numeric(x)[4:6]
  
sorted_values <- sort(values_4_to_6, decreasing = TRUE)
third_largest <- sorted_values[3]
if (third_largest< 0.3) {
    return(NA)}else{

index<-which(values_4_to_6==third_largest)[1]
return(names_1_to_3[index])}

})
combined_cell_propo<-as.data.frame(combined_cell_propo)

combined_cell_propo$top1<-top1
combined_cell_propo$top2<-top2
combined_cell_propo$top3<-top3



combined_cell_propo$top1<-gsub("B/plasma_cells","B.plasma_cells",combined_cell_propo$top1)
combined_cell_propo$top2<-gsub("B/plasma_cells","B.plasma_cells",combined_cell_propo$top2)
combined_cell_propo$top3<-gsub("B/plasma_cells","B.plasma_cells",combined_cell_propo$top3)

   
rm_combined1<-which(combined_cell_propo$CellType1 %in% c("Epithelial_cells","Erythrocytes","ILC","Mesothelial","Myofibroblasts","Pericytes"))
rm_combined2<-which(combined_cell_propo$CellType2 %in% c("Epithelial_cells","Erythrocytes","ILC","Mesothelial","Myofibroblasts","Pericytes"))

rm_combined<-c(rm_combined1,rm_combined2)
combined_cell_propo<-combined_cell_propo[-rm_combined,]
write.table(combined_cell_propo,"./simulation/pan_cancer/combined_cell_propo_summary.txt",sep="\t")


sam_acc_res<-data.frame()
ct<-0.15
top3<-0.15
top=1
X2_2<-X2
X2_2_m <- apply(X2_2, 1, max) 
rms_v1<-which(X2_2_m <ct)
rms_v2<-which(rowSums(X2_2>0.3)>2)
rms_v<-c(rms_v1,rms_v2)
print(table((rowSums(X2_2>0.2))))
X3<-X2_2[!rownames(X2_2)%in% names(rms_v),]
ind_matrix_v3<-combined_cell_propo[!rownames(combined_cell_propo)%in% names(rms_v),]

#####accuracy for overall
hit_res<-vector()
hit_names_res<-vector()
all_max<-vector()
for (ih in 1:nrow(X3)){
tmp_p<-as.numeric(X3[ih,])
names(tmp_p)<-colnames(X3)
sorted_tmp_p <- sort(tmp_p, decreasing = TRUE)
nm<-names(sorted_tmp_p)[which(sorted_tmp_p==sorted_tmp_p[top])]
tmp_r<-ind_matrix_v3[which(rownames(ind_matrix_v3)==rownames(X3)[ih]),]
hit<-ifelse(nm%in% c(tmp_r$top1,tmp_r$top2,tmp_r$top3),1,0)
names(hit)<-rownames(X3)[ih]
hit_names<-ifelse(nm%in% c(tmp_r$top1,tmp_r$top2,tmp_r$top3),nm,NA)
hit_names_res<-append(hit_names_res,hit_names)
hit_res<-append(hit_res,hit)
all_max<-append(all_max,nm)
}
hit<-names(which(hit_res==1))
miss<-names(which(hit_res==0))
hit_test<-ind_matrix_v3[hit,]
miss_test<-ind_matrix_v3[miss,]

hit_test$cn<-apply(hit_test[,7:9],1,function(x) length(which(!is.na(x))))
miss_test$cn<-apply(miss_test[,7:9],1,function(x) length(which(!is.na(x))))


table(hit_test$cn)
table(miss_test$cn)


#####accuracy for individual dataset
X3<-data.frame(X3)
X3$sample_info<-sapply(strsplit(rownames(X3), "_"), function(x) paste(x[-length(x)], collapse = "_"))
ind_matrix_v3$sample_info<-sapply(strsplit(rownames(ind_matrix_v3), "_"), function(x) paste(x[-length(x)], collapse = "_"))
#ind_matrix_v3<-ind_matrix_v[!rownames(ind_matrix_v)%in% names(rms_v),]
print(dim(X3))
acc_res<-data.frame()
n=1
for (sam in (unique(X3$sample_info))){
print(sam)
X3_tmp<-X3%>%filter(sample_info==sam)
ind_matrix_v3_tmp<-ind_matrix_v3%>%filter(sample_info==sam)
ind_matrix_v3_tmp<-ind_matrix_v3_tmp[rownames(X3_tmp),]
hit_res<-vector()
hit_names_res<-vector()
all_max<-vector()
for (ih in 1:nrow(X3_tmp)){
tmp_p<-as.numeric(X3_tmp[ih,-ncol(X3_tmp)])
names(tmp_p)<-colnames(X3_tmp)[-ncol(X3_tmp)]
sorted_tmp_p <- sort(tmp_p, decreasing = TRUE)

nm<-names(sorted_tmp_p)[which(sorted_tmp_p==sorted_tmp_p[1])]
tmp_r<-ind_matrix_v3_tmp[which(rownames(ind_matrix_v3_tmp)==rownames(X3_tmp)[ih]),]
hit<-ifelse(nm%in% c(tmp_r$top1,tmp_r$top2,tmp_r$top3),1,0)
names(hit)<-rownames(X3_tmp)[ih]
hit_names<-ifelse(nm%in% c(tmp_r$top1,tmp_r$top2,tmp_r$top3),nm,NA)
hit_names_res<-append(hit_names_res,hit_names)
hit_res<-append(hit_res,hit)
all_max<-append(all_max,nm)
}
hit<-names(which(hit_res==1))
miss<-names(which(hit_res==0))
hit_test<-ind_matrix_v3_tmp[hit,]
miss_test<-ind_matrix_v3_tmp[miss,]
hit_test$cn<-apply(hit_test[,4:6],1,function(x) length(which(as.numeric(x)>0)))
miss_test$cn<-apply(miss_test[,4:6],1,function(x) length(which(as.numeric(x)>0)))
print(table(hit_test$cn))
print(table(miss_test$cn))
if(length(table(miss_test$cn))==3 & length(table(hit_test$cn))==3){
acc_res[n,1]<-sam
acc_res[n,2]<-matrix(table(hit_test$cn)/(table(hit_test$cn)+table(miss_test$cn)))[1,1]
acc_res[n,3]<-matrix(table(hit_test$cn)/(table(hit_test$cn)+table(miss_test$cn)))[2,1]
acc_res[n,4]<-matrix(table(hit_test$cn)/(table(hit_test$cn)+table(miss_test$cn)))[3,1]
acc_res[n,5]<-dim(X3_tmp)[1]
n=n+1
}}


overall_acc<-table(hit_res)[2]/length(all_max)
print(paste0("overall accuracy: ",overall_acc))
hit_res<-data.frame(hit_res)

hit_res$sample_info<-sapply(strsplit(rownames(hit_res), "_"),function(x) paste(x[-length(x)], collapse = "_"))
sam_percent<-data.frame(tapply(hit_res[,1],as.factor(hit_res$sample_info),sum)/tapply(hit_res[,1],as.factor(hit_res$sample_info),length))
colnames(sam_percent)<-"acc"
sam_percent$condi<-paste0(cv_i,"_",ct,"_",top)
sam_acc_res<-rbind(sam_acc_res,sam_percent)
missing_names<-setdiff(names(table(all_max)),names(table(hit_names_res)))
hit_names_res<-append(hit_names_res,missing_names)
precision_s1<-table(hit_names_res)/table(all_max)

print("########precision_s1 for validation###########")
print(data.frame(precision_s1))


####calculate rmse######
X2_2<-X2
X2_2_m <- apply(X2_2, 1, max) 
rms_v1<-which(X2_2_m <0.15)
rms_v2<-which(rowSums(X2_2>0.3)>2)
rms_v<-c(rms_v1,rms_v2)
print(table((rowSums(X2_2>0.2))))
X3<-X2_2[!rownames(X2_2)%in% names(rms_v),]
colnames(X3)<-gsub("B.plasma_cells","B/plasma_cells",colnames(X3))
cell_type_matrix<-read.table("./simulation/pan_cancer/pan_sample_all_cell_type_summary.txt",sep="\t",header=T,check.names=F)

ind_matrix_v3<-cell_type_matrix[!rownames(cell_type_matrix)%in% names(rms_v),
intersect(colnames(X3),colnames(cell_type_matrix))]

X3<-data.frame(X3[,intersect(colnames(X3),colnames(cell_type_matrix))])
X3$sample_info<-sapply(strsplit(rownames(X3), "_"), function(x) paste(x[-length(x)], collapse = "_"))
ind_matrix_v3$sample_info<-sapply(strsplit(rownames(ind_matrix_v3), "_"), function(x) paste(x[-length(x)], collapse = "_"))
colnames(X3)<-gsub("B.plasma_cells","B/plasma_cells",colnames(X3))
rmse_res<-data.frame()
rmse_res_cell<-data.frame()
n=1
m=1
for (sam in (unique(X3$sample_info))){
print(sam)
X3_tmp<-X3%>%filter(sample_info==sam)
ind_matrix_v3_tmp<-ind_matrix_v3%>%filter(sample_info==sam)
row_in<-intersect(rownames(X3_tmp),rownames(ind_matrix_v3_tmp))
ind_matrix_v3_tmp<-ind_matrix_v3_tmp[row_in,]
X3_tmp<-X3_tmp[row_in,]
print(identical(rownames(X3_tmp),rownames(ind_matrix_v3_tmp)))
print(identical(colnames(X3_tmp),colnames(ind_matrix_v3_tmp)[1:10]))
sum_residuals<-0
for (i in 1:(ncol(X3_tmp)-1)){
residuals <- X3_tmp[,i] - ind_matrix_v3_tmp[,i]

# Square the residuals
squared_residuals <- residuals^2
mean_squared_residuals <- mean(squared_residuals)
rmse_res_cell[m,1]<-sam
rmse_res_cell[m,2]<-colnames(X3_tmp)[i]
rmse_res_cell[m,3]<-sqrt(mean_squared_residuals)
sum_residuals<-sum_residuals+mean_squared_residuals
m=m+1
}
# Take the square root to get RMSE
#rmse <- sqrt(sum_residuals)/ncol(X3_tmp)
rmse <- sqrt(sum_residuals/(ncol(X3_tmp)-1))
rmse_res[n,1]<-sam
rmse_res[n,2]<-rmse
n=n+1}
rmse_res_cell<-rmse_res_cell%>%filter(!V1 %in% c("RCC1","RCC3","C84"))
rmse_average <- aggregate(V3 ~ V2, data = rmse_res_cell, FUN = mean)

###########calculate pcc###################################
X3<-data.frame(X2_2[!rownames(X2_2)%in% names(rms_v),])
colnames(X3)<-gsub("B.plasma_cells","B/plasma_cells",colnames(X3))
ind_matrix_v3<-cell_type_matrix[rownames(X3),]
columns_to_add1 <- setdiff(union(colnames(X3), colnames(ind_matrix_v3)),colnames(X3))
X3[, columns_to_add1] <- 0

X3$sample_info<-sapply(strsplit(rownames(X3), "_"), function(x) paste(x[-length(x)], collapse = "_"))
ind_matrix_v3$sample_info<-sapply(strsplit(rownames(ind_matrix_v3), "_"), function(x) paste(x[-length(x)], collapse = "_"))

ind_matrix_v3<-ind_matrix_v3[,colnames(X3)]
ind_matrix_v4<-combined_cell_propo[rownames(ind_matrix_v3),]
ind_matrix_v4$cn<-apply(ind_matrix_v4[,4:6],1,function(x) length(which(as.numeric(x)>0)))
ind_matrix_v4$sample_info<-sapply(strsplit(rownames(ind_matrix_v4), "_"), function(x) paste(x[-length(x)], collapse = "_"))
ind_matrix_v4$cell.names<-rownames(ind_matrix_v4)
pcc_res<-data.frame()
n=1
pcc_average_all<-data.frame()
for (sam in (unique(X3$sample_info))){
print(sam)
X3_tmp<-X3%>%filter(sample_info==sam)
print(dim(X3_tmp))
ind_matrix_v4_tmp<-ind_matrix_v4%>%filter(sample_info==sam)
ind_matrix_v3_tmp<-ind_matrix_v3%>%filter(sample_info==sam)
ind_matrix_v3_tmp<-ind_matrix_v3_tmp[rownames(X3_tmp),]

pcc_all<-vector()
for (ri in 1:nrow(X3_tmp)){
pcc<-cor(as.numeric(X3_tmp[ri,-ncol(X3_tmp)]) ,as.numeric(ind_matrix_v3_tmp[ri,-ncol(X3_tmp)]))
pcc_all<-append(pcc_all,pcc)
}
pcc_all<-data.frame(pcc_all)
pcc_all$cell.names<-rownames(X3_tmp)
pcc_all<-merge(pcc_all,ind_matrix_v4_tmp[,10:12],by="cell.names")
pcc_average<- aggregate(pcc_all ~ cn, data = pcc_all, FUN = mean)
pcc_average$sample_info<-sam
pcc_average_all<-rbind(pcc_average_all,pcc_average)
pcc_res[n,1]<-sam
pcc_res[n,2]<-mean(pcc_all$pcc_all)
n=n+1
}

########calculate auc###############
X2_m <- apply(X2, 1, max) 
rms_v1<-which(X2_m <0.15)
rms_v2<-which(rowSums(X2>0.15)>3)
rms_v<-c(rms_v1,rms_v2)
X2_s1<-data.frame(X2[!rownames(X2)%in% names(rms_v),])
X2_s1 <- t(apply(X2_s1, 1, function(row) {
  sorted_row <- sort(row, decreasing = TRUE)
  row[row %in% sorted_row[1] & sorted_row[1] > 0.15] <- 1
  row[row %in% sorted_row[2] & sorted_row[2] > 0.15] <- 1
  row[row %in% sorted_row[3] & sorted_row[3] > 0.15] <- 1
  row[row != 1] <- 0
  return(row)
}))
ind_matrix_v_s1<-ind_matrix_v[!rownames(ind_matrix_v)%in% names(rms_v),]
auc_res_s1_tmp<-data.frame()
print(dim(X2_s1))
for (xc in 1: ncol(X2_s1)){
if(length(levels(as.factor(X2_s1[,xc])))==2 & length(levels(as.factor(ind_matrix_v_s1[,xc])))==2){
m_X_s1<-confusionMatrix(as.factor(X2_s1[,xc]), as.factor(ind_matrix_v_s1[,xc]),positive="1")  
auc_res_s1_tmp[xc,1]<-colnames(X2_s1)[xc]
auc_res_s1_tmp[xc,2]<- m_X_s1$byClass["F1"]
auc_res_s1_tmp[xc,3]<-  mltools::mcc(X2_s1[,xc], ind_matrix_v_s1[,xc])
#auc_res_tmp[xc,4]<- auc_v
auc_res_s1_tmp[xc,4]<- m_X_s1$byClass["Precision"]
}else{
auc_res_s1_tmp[xc,1]<-colnames(X3_s1)[xc]
#value_tmp<-paste0(unique(X3[,xc]),";",unique(ind_matrix_v3[,xc]))
auc_res_s1_tmp[xc,2]<-NA
auc_res_s1_tmp[xc,3]<- mltools::mcc(X2_s1[,xc], ind_matrix_v_s1[,xc])
auc_res_s1_tmp[xc,4]<- NA

}
auc_res_s1_tmp[xc,5]<-NA
auc_res_s1_tmp[xc,6]<-paste0(cv_i)
auc_res_s1_tmp[xc,7]<-length(which(ind_matrix_v_s1[,xc]==1))
}

##########public spatial validation data###############
file_folder1<-list.files("./spatial_public_data/GSE203612_RAW/st_data/")
file_folder2<-list.files("./spatial_public_data/34493872_BC/raw_count_matrices/raw_count_matrices")
file_folder3<-list.files("./spatial_public_data/GSE144239_cSCC/st_data/")
file_folder4<-list.files("./spatial_public_data/GSE190811_BC/st_data/")
file_folder5<-list.files("./spatial_public_data/GSE175540_RC/st_data/")
file_folder6<-list.files("./spatial_public_data/PMC8683021_LiC/st_data/")
file_name5<-file_folder5
file_name6<-file_folder6
file_name4<-file_folder4[3:4]
file_name3<-file_folder3

file_name2<-str_split_fixed(file_folder2, pattern="_",n=2)
file_name2<-file_name2[,1]
file_name1<-file_folder1
file_name<-c(file_name1,file_name2,file_name3,file_name4,file_name5,file_name6)
cnv_res<-data.frame()
for (i in 1:length(file_name)){
print(file_name[i])
cnv_data<-read.table(paste0("./simulation/public_data/cnv/",file_name[i],"_0803_copykat_prediction.txt"),header=T,sep="\t")
cnv_data<-cnv_data%>%filter(copykat.pred!="not.defined"& grepl(paste0("^",file_name[i]), cell.names))
cnv_res<-rbind(cnv_res,cnv_data)
}
write.table(cnv_res,"./simulation/public_data/sample_cnvres_0803.txt",sep="\t")


cnv_res<-read.table("./simulation/public_data/sample_cnvres_0803.txt",sep="\t",header=T)
file_folder<-list.files("./spatial_public_data/gsva_res/")
file_name<-file_folder
gsva_res<-data.frame()
for (i in 1:length(file_name)){
print(file_name[i])
gsva_data<-read.table(paste0("./spatial_public_data/gsva_res/",file_name[i]),header=T,sep="\t")
print(dim(gsva_data))
#print(colnames(gsva_data))
gsva_res<-rbind(gsva_res,gsva_data)
}

gsva_res$cell.names<-rownames(gsva_res)
cnv_res$cell.names<-gsub("CID4290A","CID4290",cnv_res$cell.names)
com_res<-merge(cnv_res,gsva_res,by="cell.names")
com_res$sample_info<-sapply(strsplit(com_res$cell.names, "_"), function(x) paste(x[-length(x)], collapse = "_"))
data_fea<-com_res[,c(1:54,69)]
# write.table(com_res,"./spatial_public_data/cnv2/19_sample_cnv_gsvares.txt",sep="\t")
# com_res<-read.table("./spatial_public_data/cnv2/19_sample_cnv_gsvares.txt",sep="\t",header=T)
# data_fea<-all_data[,c(1:54,69)]
data_fea$copykat.pred2<-ifelse(data_fea$copykat.pred=="aneuploid",1,0)
data_fea[is.na(data_fea)]<-0
data_fea2<-data.frame()
for (sam in unique(data_fea$sample_info)){
tmp<-data_fea%>%filter(sample_info==sam)
tmp[is.na(tmp)]<-0
tmp[,c(3:54)]<-apply(tmp[,c(3:54)],c(1,2),as.numeric)
tmp[,c(3:54)]<-apply(tmp[,c(3:54)],2,function(x) (x-min(x))/(max(x)-min(x)))
tmp[is.na(tmp)]<-0
data_fea2<-rbind(data_fea2,tmp)}
# data_fea2$class<-all_data$class
X<-data.frame()
for (sam in unique(data_fea2$sample_info)){
data_fea3<-data_fea2%>%filter(data_fea2$sample_info==sam)
Y_mat_all<-data_fea3[,c(1,3:54,56,55)]
rownames(Y_mat_all)<-Y_mat_all$cell.names
Y_mat_all<-Y_mat_all[,-1]
Y_mat_all_v1<-Y_mat_all
#Y_mat_all_v1<-Y_mat_all_v1[,-61]
Y_mat_all_v1<-Y_mat_all_v1[,-54]
Y_mat_all_v<-as.matrix(Y_mat_all_v1)
Y_mat_all_v<-apply(Y_mat_all_v,c(1,2),as.numeric)
Y<-Y_mat_all_v
set.seed(123123)
X_tmp<-Y%*%t(A)%*%solve((A %*% t(A)) +lambda * diag(nrow(A)))
X<-as.data.frame(X)
X[X<0]=0
X<-rbind(X,X_tmp)
}

X2=X
X2<-t(apply(X2,1,function(x) x/sum(x)))
X2[is.nan(X2)] <- 0
X2_m <- apply(X2, 1, max) 
####cutoff1######
rms_v1<-which(X2_m <0.15)
####cutoff2######
rms_v2<-which(rowSums(X2>0.3)>2)
rms_v<-c(rms_v1,rms_v2)
X2_s1<-data.frame(X2[!rownames(X2)%in% names(rms_v),])

X2_top<-data.frame(X2[!rownames(X2)%in% names(rms_v),])
X2_top$top1<-unlist(apply(X2_top[,1:9],1,function(x) names(which.max(x))[1]))

for (ih in 1:nrow(X2_top)){
tmp_p<-as.matrix(X2_top[,1:9])[ih,]
sorted_tmp_p <- sort(tmp_p, decreasing = TRUE)
if(sorted_tmp_p[2]*0.83<0.15){
X2_top$top2[ih]<-NA
}else{
X2_top$top2[ih]<-names(sorted_tmp_p)[which(sorted_tmp_p==sorted_tmp_p[2])]
}}
for (ih in 1:nrow(X2_top)){
tmp_p<-as.matrix(X2_top[,1:9])[ih,]
sorted_tmp_p <- sort(tmp_p, decreasing = TRUE)
if(sorted_tmp_p[3]*0.75<0.15){
X2_top$top3[ih]<-NA
}else{
X2_top$top3[ih]<-names(sorted_tmp_p)[which(sorted_tmp_p==sorted_tmp_p[3])]
}}

X2_top$spot_class<-ifelse("Malignant" %in% X2_top[,c("top1","top2","top3")],"tumor","non_tumor")

X2_top$spot_class<-apply(X2_top[,c("top1","top2","top3")],1,
function(x) ifelse("Malignant" %in% x & length(which(is.na(x))=="TRUE")==2,"tumor","tme"))
X2_top$spot_class<-apply(X2_top[,c("top1","top2","top3","spot_class")],1,
function(x) ifelse("Malignant" %in% x & x[4]=="tme","tumor_tme",x[4]))
X2_top$sample_info<-sapply(strsplit(rownames(X2_top), "_"), function(x) paste(x[-length(x)], collapse = "_"))


# identical(rownames(X2_top),rownames(X2_s2))
# X2_all<-cbind(X2_top,X2_s2)
X2_all<-X2_top
X2_all_res<-X2_all[,c("top1","top2","top3","spot_class","sample_info","Malignant")]
X2_uniden<-data.frame(X2[rownames(X2)%in% names(rms_v),])
X2_uniden$spot_class<-X2_uniden$top3<-X2_uniden$top2<-X2_uniden$top1<-"unidentified"
X2_uniden$sample_info<-sapply(strsplit(rownames(X2_uniden), "_"), function(x) paste(x[-length(x)], collapse = "_"))
X2_uniden_res<-X2_uniden[,c("top1","top2","top3","spot_class","sample_info","Malignant")]

X2_res<-rbind(X2_all_res,X2_uniden_res)
X2_res$Malignant[is.na(X2_res$Malignant)]<-0
table(X2_res$sample_info,X2_res$top1)
table(X2_res$sample_info,X2_res$spot_class)
X2_res<-read.table("./spatial_public_data/publicdata_pre_res_all.txt",sep="\t",header=T)
X2_ori_res<-read.table("./spatial_public_data/publicdata_pre_original_res_all.txt",sep="\t",header=T)

