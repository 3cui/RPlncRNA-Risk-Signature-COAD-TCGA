exp=read.csv("TCGA-COAD_mrna_expr_counts.csv",row.names = 1)
exp1=exp[rownames(exp)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]
nrow(exp)
AA=as.data.frame(rownames(exp))
exp1 = exp[rowSums(exp)>0,]
nrow(exp1)
AA=as.data.frame(rownames(exp1))
rownames(AA)=AA$`rownames(exp1)`
AA1=AA[rownames(AA)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]
table(substr(colnames(exp1),14,16))
Tumor <- grep('01A',colnames(exp1))
Tumor  
Tumor_exp1=exp1[,Tumor]
Normal <- grep('11A',colnames(exp1))
Normal 
Normal_exp1=exp1[,Normal]
counts_exp1=cbind(Normal_exp1,Tumor_exp1)
table(substr(colnames(counts_exp1),14,16))
Tumor <- grep('01A',colnames(counts_exp1))
Tumor
Group <- factor(c(rep("Normal",times=41),rep("Tumor",times=465)))

Data <- data.frame(row.names = colnames(counts_exp1), 
                   group = Group)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_exp1,
                              colData = Data,
                              design = ~ group)
dds2 <- DESeq(dds)
res <- results(dds2, contrast=c("group", "Tumor", "Normal"))#肿瘤在前，对照在后
res <- res[order(res$pvalue),]
summary(res)
my_result <- as.data.frame(res)
my_result <- na.omit(my_result)
library(dplyr)
my_result$Gene_symbol<-rownames(my_result)
my_result <- my_result %>% dplyr::select('Gene_symbol',colnames(my_result)[1:dim(my_result)[2]-1],everything())
rownames(my_result) <- NULL

my_result$regulate <- ifelse(my_result$padj > 0.05, "unchanged",
                             ifelse(my_result$log2FoldChange > log2(1.5), "up-regulated",
                                    ifelse(my_result$log2FoldChange < -log2(1.5), "down-regulated", "unchanged")))
table(my_result$regulate)
rownames(my_result)=my_result$Gene_symbol
write.csv(my_result,file= "1.2DEG_deseq2(mRNAlog2(1.5)).csv")

my_result2=my_result[rownames(my_result)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]
my_result3=my_result2[my_result2$regulate!="unchanged",]


exp=read.csv("TCGA-COAD_lncrna_expr_counts.csv",row.names = 1)
nrow(exp)
exp1 = exp[rowSums(exp)>0,]
nrow(exp1)
exp1 = exp1[apply(exp1, 1, function(x) sum(x > 0) > 0.5*ncol(exp1)), ]
nrow(exp1)
table(substr(colnames(exp1),14,16))
Tumor <- grep('01A',colnames(exp1))
Tumor  
Tumor_exp1=exp1[,Tumor]
Normal <- grep('11A',colnames(exp1))
Normal
Normal_exp1=exp1[,Normal]
counts_exp1=cbind(Normal_exp1,Tumor_exp1)
table(substr(colnames(counts_exp1),14,16))
Tumor <- grep('01A',colnames(counts_exp1))
Tumor
Group <- factor(c(rep("Normal",times=41),rep("Tumor",times=465)))
Data <- data.frame(row.names = colnames(counts_exp1),
                   group = Group)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_exp1,
                              colData = Data,
                              design = ~ group)
dds2 <- DESeq(dds)
res <- results(dds2, contrast=c("group", "Tumor", "Normal"))
res <- res[order(res$pvalue),]
summary(res)
my_result <- as.data.frame(res)
my_result <- na.omit(my_result)
library(dplyr)
my_result$Gene_symbol<-rownames(my_result)
my_result <- my_result %>% dplyr::select('Gene_symbol',colnames(my_result)[1:dim(my_result)[2]-1],everything())
rownames(my_result) <- NULL

my_result=read.csv("1.1DEG_deseq2(lncRNA1).csv",row.names = 1)
my_result$regulate <- ifelse(my_result$padj > 0.05, "unchanged",
                             ifelse(my_result$log2FoldChange > log2(1.5), "up-regulated",
                                    ifelse(my_result$log2FoldChange < -log2(1.5), "down-regulated", "unchanged")))
rownames(my_result)=my_result$Gene_symbol
my_result=read.csv("1.1DEG_deseq2(lncRNA1).csv",row.names = 1)
table(my_result$regulate)
my_result=my_result[my_result$regulate!="unchanged",]####只要有表达的基因
exp=read.csv("TCGA-COAD_lncrna_expr_fpkm.csv",row.names = 1)
table(substr(colnames(exp),14,16))
Tumor <- grep('01A',colnames(exp))
Tumor_exp=exp[,Tumor]
Normal <- grep('11A',colnames(exp))
Normal 
Normal_exp=exp[,Normal]
fpkm_exp=cbind(Normal_exp,Tumor_exp)
exp=read.csv("TCGA-COAD_mrna_expr_fpkm.csv",row.names = 1)###泛凋亡基因不只是lncRNA，所以要用全部的fpkm数据
table(substr(colnames(exp),14,16))
Tumor <- grep('01A',colnames(exp))
Tumor  
Tumor_exp=exp[,Tumor]
Normal <- grep('11A',colnames(exp))
Normal 
Normal_exp=exp[,Normal]
fpkm_exp=cbind(Normal_exp,Tumor_exp)
gene=read.csv("1.2DEG_deseq2(mRNAlog2(1.5)).csv",row.names = 1)
gene1=gene[rownames(gene)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]
fpkm_exp1=fpkm_exp[rownames(gene1),]####用全部的凋亡基因来做
write.csv(fpkm_exp1,"3.2fpkm表达矩阵中提取出泛凋亡基因.csv")

exp1=read.csv("3.1fpkm表达矩阵中提取出差异基因.csv",row.names = 1)
exp2=read.csv("3.2fpkm表达矩阵中提取出泛凋亡基因.csv",row.names = 1)
exp3=t(exp1)
exp3=as.data.frame(exp3)
range(exp3)
exp4=t(exp2)
exp4=as.data.frame(exp4)
range(exp4)
exp5=cbind(exp3,exp4)
range(exp5)
data_use=exp5
colnames(exp4)#就是用的32个凋亡基因
target_gene <- c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1")
result <- c()
gene_list <- colnames(data_use)[1:3546]#进行相关性分析的范围【我们是选择全部的基因】
for (gene in gene_list) {
  for(target in target_gene){
    message(gene,target)#打印【报出信息】
    cor_R <- cor.test(x = data_use[,target],y = data_use[,gene],method = 'pearson')###基因近似正态分布
    P.value<-cor_R[["p.value"]]
    cor=cor_R[["estimate"]][["cor"]]
    temp_result<-c(gene,target,cor,P.value)#让这个变量一直在变
    result <- rbind(result, temp_result)
  }
  
}
colnames(result)=c("lncRNA","mRNA","cor_R","cor_P")
result=as.data.frame(result)
result$cor_R=as.numeric(result$cor_R)
result$cor_P=as.numeric(result$cor_P)
table(result$cor_R > 0.4)
table(result$cor_P < 0.05)
result=read.csv("4.1 32个凋亡基因的分析.csv")
index <- result$cor_P < 0.05
index1<-result$cor_R>0.4
result_filt <- result[index,]
result_filt=result_filt[result_filt$cor_R>0.4|result_filt$cor_R< -0.4,]
library(dplyr)
result_filt1 <- result_filt[!duplicated(result_filt$lncRNA),]
rownames(result_filt1)=result_filt1$lncRNA
result_filt1=result_filt1[,-1]

exp=read.csv("output_mRNA_lncRNA_expr/TCGA-COAD_lncrna_expr_fpkm.csv",row.names = 1)
table(substr(colnames(exp),14,16))
Tumor <- grep('01A',colnames(exp))
Tumor_exp=exp[,Tumor]
Normal <- grep('11A',colnames(exp))
Normal 
Normal_exp=exp[,Normal]
fpkm_exp=cbind(Normal_exp,Tumor_exp)
fpkm_exp1=fpkm_exp[rownames(result_filt1),]
range(fpkm_exp1)
fpkm_exp1=log2(fpkm_exp1+1)
able(substr(colnames(fpkm_exp1),14,16))
Tumor <- grep('01A',colnames(fpkm_exp1))
Tumor_lnc=fpkm_exp1[,Tumor]
library(stringr)
expr_use = Tumor_lnc[,sort(colnames(Tumor_lnc))]
k = !duplicated(str_sub(colnames(expr_use),1,12))
table(k)
expr_use = expr_use[,k] 
colnames(expr_use) <- substr(colnames(expr_use),1,12)
survival_data=survival_data[c(1,5,2,3,4,6,7,8,9)]
table(substr(rownames(survival_data),14,16))
Tumor <- grep('01A',rownames(survival_data))
survival_data=survival_data[Tumor,]
survival_data = survival_data[sort(rownames(survival_data)),]
library(stringr)
k = !duplicated(str_sub(rownames(survival_data),1,12))
rownames(survival_data) <- substr(rownames(survival_data),1,12)
library(stringr)
rownames(survival_data)=str_replace_all(rownames(survival_data),"-",".")
datExpr <- t(expr_use)
datExpr=as.data.frame(datExpr)
dim(datExpr)
dim(survival_data)
head(rownames(datExpr))
head(rownames(survival_data))
s = intersect(rownames(datExpr),rownames(survival_data));length(s)
datExpr1 = datExpr[s,]
survival_data2 = survival_data[s,]
dim(datExpr1)
dim(survival_data2)
identical(rownames(survival_data2),rownames(datExpr1))
datExpr <- cbind(survival_data2,datExpr1)
lnc=read.csv("4.2去重复后的数据框(0.3).csv",row.names = 1)
lnc=lnc[lnc$cor_R>0.4,]
library(stringr)
lnc1=str_replace(rownames(lnc),"-",".")
lnc1=as.data.frame(lnc1)
datExpr1=datExpr[c(1:2)]
K=intersect(lnc1$lnc1,colnames(datExpr))
datExpr2=datExpr[K]
datExpr=cbind(datExpr1,datExpr2)
library(caret)
set.seed(100)
sam<- createDataPartition(datExpr$status, p = .7,list = FALSE)
train <- datExpr[sam,]
test <- datExpr[-sam,]
write.csv(train,"6.2训练集.csv")
write.csv(test,"6.3验证集.csv")
library(tidyverse)
library(survival)
library(rms)
library(survminer)
library(forestmodel)
library(forestplot)
library(ggplot2)
datExpr=read.csv("6.2训练集.csv",row.names = 1)
summary(datExpr$time)
table(datExpr$status)
dd <- datadist(datExpr)
options(datadist = "dd")
summary(datExpr)
surv_object <- with(datExpr, Surv(time, status))
result <- data.frame("Gene" = character(),
                     "Hazard Ratio" = numeric(),
                     "95%CI" = character(),
                     "P value" = numeric())

for (i in 3:length(colnames(datExpr))) {
  print(i)
  gene <- colnames(datExpr)[i]
  model <- coxph(surv_object ~ datExpr[[gene]], data = datExpr)
  data_use <- summary(model)
  HR <- round(data_use$coefficients[,2],2)
  CI <- paste0(round(data_use$conf.int[,3:4],2),collapse = '-')
  P_value <- round(data_use$coefficients[,5],3)
  temp_result <- data.frame("Gene" = colnames(datExpr)[i],
                            "Hazard Ratio" = HR,
                            "95%CI" = CI,
                            "P value" = P_value)
  result <- rbind(result, temp_result)
}
result_use <- result
result_use$significant <- ifelse(result_use$P.value < 0.05,"significant","Not significant")
table(result_use$significant)
colnames(result_use)
result_use1=result_use[result_use$P.value < 0.05,]
result_use1=result_use1[result_use1$Hazard.Ratio!="0",]
survival_data=datExpr[,c(1,2)]
datExpr1 <- datExpr[,-c(1:2)]
datExpr1=t(datExpr1)
datExpr1=as.data.frame(datExpr1)
datExpr2=datExpr1[result_use1$Gene,]
datExpr2=t(datExpr2)
datExpr2=as.data.frame(datExpr2)
datExpr=cbind(survival_data,datExpr2)
write.csv(datExpr,"7.3单因素cox结果与临床和表达矩阵的结合.csv")
library(glmnet)
library(corrplot)
library(car)
library(survival)
library(dplyr)
summary(datExpr$time)
table(datExpr$status)
min(datExpr$time)
datExpr <- datExpr %>% filter(time > 0)
data_use <- datExpr
nrow(data_use)
ncol(data_use)
time <- data_use[,"time"]
status <- data_use[,"status"]
X <- as.matrix(data_use[3:25])
Y <- as.matrix(Surv(time, status))
lasso <- glmnet(X,Y,family = "cox",alpha = 1)
set.seed(1000)
lasso_cv <- cv.glmnet(X,Y,family = "cox",alpha = 1,nfolds = 5)
plot(lasso_cv) 
lambda <- lasso_cv$lambda.min
coef_lasso_cv <- coef(lasso, s = lambda)
coef_lasso_cv
exp(coef_lasso_cv)
lasso_lncRNA=read.csv("8.4lasso.csv",row.names = 1)
lasso_lncRNA$gene=rownames(lasso_lncRNA)
lasso_lncRNA1=lasso_lncRNA[lasso_lncRNA$X1!="0",]
bind=read.csv("7.3单因素cox结果与临床和表达矩阵的结合.csv",row.names = 1)
bind=datExpr
bind1=bind[,c(1,2)]
bind2=bind[,rownames(lasso_lncRNA1)]
bind3=cbind(bind1,bind2)
datExpr=read.csv("8.5lasso回归结果与临床和表达矩阵的结合.csv",row.names = 1)
colnames(datExpr)
datExpr1=datExpr[,-c(1,2)]
formula <- as.formula(paste("Surv(time, status) ~", paste(colnames(datExpr1), collapse = " + ")))
multi_cox_model <- coxph(formula, data = datExpr)
data_use <- summary(multi_cox_model)
ulti_cox_HR <- round(data_use$coefficients[,2],2)
multi_cox_CI2.5 <- round(data_use$conf.int[,3],2)
multi_cox_CI97.5 <- mul_CI95<-round(data_use$conf.int[,4],2)
multi_cox_CI <- paste0('(',multi_cox_CI2.5,'-',multi_cox_CI97.5,')')
multi_cox_P_value <- round(data_use$coefficients[,5],3)
Variable <- row.names(data.frame(data_use$coefficients))
multi_cox_result<- data.frame(Variable,multi_cox_HR,multi_cox_CI2.5,multi_cox_CI97.5,multi_cox_CI,multi_cox_P_value)
formula <- as.formula(paste("Surv(time, status) ~", paste(rownames(multi_cox_result1), collapse = " + ")))
multi_cox_model <- coxph(formula, data = datExpr)
summary(multi_cox_model)
data_use <- summary(multi_cox_model)
multi_cox_HR <- round(data_use$coefficients[,2],2)
multi_cox_CI2.5 <- round(data_use$conf.int[,3],2)
multi_cox_CI97.5 <- mul_CI95<-round(data_use$conf.int[,4],2)
multi_cox_CI <- paste0('(',multi_cox_CI2.5,'-',multi_cox_CI97.5,')')
multi_cox_P_value <- data_use$coefficients[,5]
Variable <- row.names(data.frame(data_use$coefficients))
multi_cox_result<- data.frame(Variable,multi_cox_HR,multi_cox_CI2.5,multi_cox_CI97.5,multi_cox_CI,multi_cox_P_value)
write.csv(multi_cox_result, file = "11.4多因素cox_results(基因筛选).csv")
multi_cox_model <- coxph(formula, data = datExpr)
summary(multi_cox_model)
forest_model(multi_cox_model) 
datExpr=read.csv("11.1训练集分数计算.csv",row.names = 1)
datExpr$D1=datExpr$LINC01133*datExpr$C1
datExpr$D2=datExpr$FOXD3.AS1*datExpr$C2
datExpr$D3=datExpr$AP001066.1*datExpr$C3
datExpr$D4=datExpr$AP003555.1*datExpr$C4
datExpr$risk.score=datExpr$D1+datExpr$D2+datExpr$D3+datExpr$D4
library(tidyverse)
library(survival)
library(rms)
library(survminer)
library(forestmodel)
library(forestplot)
library(ggplot2)
mydata=datExpr
head(mydata)
res.cut <- surv_cutpoint(mydata, time = "time", event = "status",variables = c("risk.score"))###基因名字
summary(res.cut)
plot(res.cut, "risk.score", palette = "npg")

res.cat <- surv_categorize(res.cut)
head(res.cat)
table(res.cat$risk.score)
write.csv(res.cat,"11.1训练集cox_risk.score__KM_group.csv")
identical(rownames(res.cat),rownames(mydata))
mydata$group_risk.score=res.cat$risk.score
write.csv(mydata,"11.2训练集cutoff分组完成与表达矩阵结合.csv")
mydata=mydata
mydata$group=ifelse(mydata$risk.score>median(mydata$risk.score),"high","low")
table(mydata$group)
write.csv(mydata,"11.3训练集中位数分组完成与表达矩阵结合.csv")
mydata=read.csv("11.3训练集中位数分组完成与表达矩阵结合.csv",row.names = 1)
library(survival)
library(survminer)
colnames(mydata)
fit <- survfit(Surv(time, status) ~group, data = mydata)
library(ggplot2)
ggsurvplot(fit, data = mydata, 
           ggtheme = theme_bw(),
           risk.table = TRUE, 
           risk.table.col = "strata",
           conf.int = FALSE,
           pval = TRUE)
library(timeROC)
library(survival)
library(tidyverse)
exp_sur=mydata
exp_sur$time <- exp_sur$time/365
ROC3 <- timeROC(T=exp_sur$time,   #结局时间
                delta=exp_sur$status,   #结局指标
                marker=exp_sur$risk.score,   #预测变量
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1, 3, 5),   #时间点，选取1年，3年和5年的生存率
                iid=TRUE)
plot(ROC3,
     time=1, col="red")   
plot(ROC3,
     time=3, col="green", add=TRUE)   
plot(ROC3,
     time=5, col="blue", add=TRUE)
legend("bottomright",
       c("AUC Year-1:76.83", "AUC Year-3:78.21", "AUC Year-5:71.40"),
       col=c("red", "green", "blue"),
       lty=1, lwd=2)   
my_result=read.csv("12.3验证集中位数分组完成与表达矩阵结合.csv",row.names = 1)
my_result1=my_result[my_result$group=="low",]
my_result2=my_result[my_result$group=="high",]
my_result3=rbind(my_result1,my_result2)#######按照风险分组重新排序!!
my_result3$patientid=1:length(my_result3$status)
my_result4=my_result3[c(1092,1089,1091)]
my_result4$risk.score=as.numeric(sort(my_result4$risk.score))####风险评分也要一定按照排序来!!，不然画图没有办法连续了，而是画了一堆点图
getwd()
write.csv(my_result4,"14.1训练集——风险三联动图矩阵1.csv")
my_result5=my_result3[c(1092,2,1)]
write.csv(my_result5,"14.2训练集——风险三联动图矩阵2.csv")
my_result5$status=ifelse(my_result5$status==0,'alive','death')
my_result5$status=factor(my_result5$status,levels = c("death","alive"))
my_result_gene=read.csv("12.3验证集中位数分组完成与表达矩阵结合.csv",row.names = 1)
my_result_gene=my_result_gene[colnames(my_result_gene)%in%c("LINC01133","FOXD3.AS1","AP003555.1","AP001066.1")]
my_result_gene=t(my_result_gene)
my_result_gene=as.data.frame(my_result_gene)
range(my_result_gene)
my_result_gene=my_result_gene[,rownames(my_result4)]
identical(rownames(my_result3),colnames(my_result_gene))
p1=ggplot(my_result4,aes(x=patientid,y=risk.score))+
  geom_point(aes(color = group))+
  scale_color_manual(values = c("#f87669","#2874C5"))+
  #geom_vline(xintercept = 107,lty = 2)+##这个是low的最后一个
  geom_vline(xintercept = 0.5*nrow(my_result5),lty = 2)+###geom_vline用于绘制垂线
  scale_x_continuous(expand=c(0.01,0))+##强制设置坐标的起始点
  labs(x = "Patients(increasing risk score)",y = "risk score")+
  theme_bw();p1
###这样连续起来是因为——我们之前按照风险评分排序了

###第二个图----
p2=ggplot(my_result5,aes(x=patientid,y=time))+
  geom_point(aes(col=status))+
  scale_color_manual(values = c("#f87669","#2874C5"))+
  #geom_vline(xintercept = 107,lty = 2)+
  geom_vline(xintercept = 0.5*nrow(my_result5),lty = 2)+
  scale_x_continuous(expand=c(0.01,0))+
  labs(x = "Patients(increasing risk score)")+
  theme_bw();p2
annotation_col = data.frame(group = my_result4$group,##分组
                            row.names = rownames(my_result4))##样本名
mycolors <- colorRampPalette(c("#2874C5","white", "#f87669"), bias = 1.2)(100)
ann_colors = list(
  group = c(low="#2874C5", high="#f87669")
)
p3=pheatmap::pheatmap(my_result_gene,#####因为annotation_col我用的是my_result4，所以样本顺序一定要保持一致
                      col= mycolors,
                      annotation_col = annotation_col,
                      annotation_colors =ann_colors,
                      scale = "row",
                      breaks = seq(-3,3,length.out = 100),###这里的参数有必要调整一下
                      show_colnames = F,
                      cluster_cols = F,
                      annotation_legend = T)
library(patchwork)
library(ggplotify)

library(rms)
library(survival)
library(foreign)
library(VIM)
library(pROC)
library(timeROC)
datExpr=read.csv("17.1 训练集临床性状cox分析矩阵.csv",row.names = 1)
colnames(datExpr)
datExpr=datExpr[c(1,2,3,4,6,10)]
colnames(datExpr)
colnames(datExpr)=c("status", "time" ,"Gender","Age","Tumor_stage","Risk_score")
aggr(datExpr, col=c('skyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(datExpr), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Percentage"))

datExpr <- na.omit(datExpr)
summary(datExpr)
summary(datExpr$time)
units(datExpr$time) <- "day"
dd <- datadist(datExpr)
options(datadist = "dd")
min(datExpr$time)
surv_object <- with(datExpr, Surv(time, status==1))
colnames(datExpr)
model <-  cph(surv_object ~ Risk_score+Age+Tumor_stage+Gender, x = TRUE ,y = TRUE, surv = TRUE, data = datExpr)
model
surv <- Survival(model)#拟合生存函数
surv_1y <- function(x)surv(365,lp=x)#一年生存函数
surv_3y <- function(x)surv(365*3,lp=x)#三年生存函数
surv_5y <- function(x)surv(365*5,lp=x)#五年生存函数
Nomogram_1 <- nomogram(model,fun = list(surv_1y,surv_3y,surv_5y),lp=F,#模型，#要放入的生存函数，#风险预测轴，T或F
                       funlabel = c('1 year survival rate','3 year survival rate','5 year survival rate'),#风险预测轴的名称
                       maxscale = 100,fun.at = c(0.1,seq(0.1,0.9,by=0.1),0.90))#分数轴总分，#风险预测轴概率取值
plot(Nomogram_1)
pred <- predict(model,datExpr,type="lp")# 使用模型预测数据
colnames(datExpr)
ROC_table <- data.frame(time = datExpr[,"time"],status = datExpr[,"status"],score = pred)
time_roc_res <- timeROC(T = ROC_table$time,
                        delta = ROC_table$status,
                        marker = ROC_table$score,
                        cause = 1,
                        weighting="marginal",
                        times = c(365, 3*365, 5*365),
                        ROC = TRUE,
                        iid = TRUE
)
time_ROC_df <- data.frame(TP_1year = time_roc_res$TP[, 1],##真阳性率
                          FP_1year = time_roc_res$FP[, 1],##假阳性率
                          TP_3year = time_roc_res$TP[, 2],
                          FP_3year = time_roc_res$FP[, 2],
                          TP_5year = time_roc_res$TP[, 3],
                          FP_5year = time_roc_res$FP[, 3]
)
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#BC1328") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",x = 0.70, y = 0.25, size = 5.5,label = paste0("AUC of 5-year survival = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#BC1328") +
  labs(x = "1-specificity", y = "Sensitivity") +
  theme(axis.text = element_text(face = "bold", size = 11, color = "black"),
        axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
        axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#0067B5") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#09891D") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#BC1328") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",x = 0.75, y = 0.20, size = 4.5,label = paste0("AUC of 1-year survival = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#0067B5") +
  annotate("text",x = 0.75, y = 0.15, size = 4.5,label = paste0("AUC of 3-year survival = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#09891D") +
  annotate("text",x = 0.75, y = 0.10, size = 4.5,label = paste0("AUC of 5-year survival = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#BC1328") +
  labs(x = "1-specificity", y = "Sensitivity") +
  theme(axis.text = element_text(face = "bold", size = 11, color = "black"),
        axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
        axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))

fit1 <- cph(Surv(time,status)~ risk.score+age+stage_copy+gender,
            time.inc = 1*365,x=T,y=T,surv = T,data=datExpr)
calibration <- calibrate(fit1, cmethod='KM', method='boot', u=1*365, m=70, B=1000)
plot(calibration,lwd=2,lty=1,errbar.col=c(rgb(0,0,0,maxColorValue=255)),xlim=c(0,1),
     ylim=c(0,1),xlab="Nomogram-Predicted Probability (%)",
     ylab="Observed OS (%)",col=c("blue"),
     subtitles = F)
fit2 <- cph(Surv(time,status)~ risk.score+gender+age+stage_copy,
            time.inc = 3*365,x=T,y=T,surv = T,data=datExpr)

calibration1 <- calibrate(fit2, cmethod='KM', method='boot', u=3*365, m=70, B=1000)
plot(calibration1,lwd=2,lty=1,errbar.col=c(rgb(0,0,0,maxColorValue=255)),xlim=c(0,1),
     ylim=c(0,1),
     col=c("green"),
     add=T)
fit3 <- cph(Surv(time,status)~ risk.score+gender+age+stage_copy,
            time.inc = 5*365,x=T,y=T,surv = T,data=datExpr)

calibration2 <- calibrate(fit3, cmethod='KM', method='boot', u=5*365, m=70, B=1000)

plot(calibration2,lwd=2,lty=1,errbar.col=c(rgb(0,0,0,maxColorValue=255)),xlim=c(0,1),
     ylim=c(0,1),
     col=c("red"),
     add=T)
legend(0.8,0.25,
       c("1-year",
         "3-year",
         "5-year"
       ),
       lty = c(1,1,1),
       lwd = c(2,2,2),
       col = c("blue","green","red"),
       bty = "0")
library(survival)
library(survminer)

setwd("D://R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
mydata <- read.csv("19.1 全部数据集临床性状cox分析矩阵.csv",row.names = 1)
group=read.csv("13.3全部中位数分组完成与表达矩阵结合.csv",row.names = 1)
identical(rownames(mydata),rownames(group))####一定要用这一步进行确认，否则后面就会出错
mydata=mydata[rownames(group),]
mydata$group=group$group
table(mydata$group)
write.csv(mydata,"23.1临床性状KM曲线的矩阵.csv")

table(mydata$stage_copy1)
mydata1=mydata[mydata$stage_copy1=="Stage I&II",]
fit <- survfit(Surv(time, status) ~group, data = mydata1)
ggsurvplot(fit, data = mydata1, 
           ggtheme = theme_bw(),
           risk.table = TRUE, 
           risk.table.col = "strata",
           conf.int = FALSE,
           pval = TRUE,
           title = "Patients with Stage I_II")

table(mydata$stage_copy1)
mydata1=mydata[mydata$stage_copy1=="Stage III&Stage IV",]
fit <- survfit(Surv(time, status) ~group, data = mydata1)
ggsurvplot(fit, data = mydata1, 
           ggtheme = theme_bw(),
           risk.table = TRUE, 
           risk.table.col = "strata",
           conf.int = FALSE,
           pval = TRUE,
           title = "Patients with Stage III_IV")


table(mydata$gender)
mydata1=mydata[mydata$gender=="male",]
fit <- survfit(Surv(time, status) ~group, data = mydata1)
ggsurvplot(fit, data = mydata1, 
           ggtheme = theme_bw(),
           risk.table = TRUE, 
           risk.table.col = "strata",
           conf.int = FALSE,
           pval = TRUE,
           title = "Patients with male")
table(mydata$gender)
mydata2=mydata[mydata$gender=="female",]
fit <- survfit(Surv(time, status) ~group, data = mydata2)
ggsurvplot(fit, data = mydata2, 
           ggtheme = theme_bw(),
           risk.table = TRUE, 
           risk.table.col = "strata",
           conf.int = FALSE,
           pval = TRUE,
           title = "Patients with female")

mydata$group_age=ifelse(mydata$age>60,"old","young")
mydata3=mydata[mydata$group_age=="old",]
fit <- survfit(Surv(time, status) ~group, data = mydata3)
ggsurvplot(fit, data = mydata3, 
           ggtheme = theme_bw(),
           risk.table = TRUE, 
           risk.table.col = "strata",
           conf.int = FALSE,
           pval = TRUE,
           title = "Patients with age >60")
mydata4=mydata[mydata$group_age=="young",]
fit <- survfit(Surv(time, status) ~group, data = mydata4)
#pdf("survival.pdf",width = 5,height = 5.5,onefile = FALSE)
library(ggplot2)
ggsurvplot(fit, data = mydata4, 
           ggtheme = theme_bw(),
           risk.table = TRUE, 
           risk.table.col = "strata",
           conf.int = FALSE,
           pval = TRUE,
           title = "Patients with age< =60")

library(ggplot2)
library(ggpubr)
mydata <- read.csv("23.1临床性状KM曲线的矩阵.csv",row.names = 1)

mydata1=mydata[mydata$group=="high",]
my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage I", "Stage IV") ,c("Stage III", "Stage II"),c("Stage IV", "Stage II"),c("Stage IV", "Stage III"))
ggboxplot(mydata1, x = "stage_copy", y = "risk.score",
          color = "stage_copy")+ 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")


library("GenomicFeatures")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#install.packages("oncoPredict")
library(oncoPredict)
#直接看 v2的版本
trainingExprData=readRDS(file='GDSC2_Expr (RMA Normalized and Log Transformed).rds')
dim(trainingExprData) #
trainingExprData[1:4,1:4]
trainingPtype=readRDS("GDSC2_Res.rds")
trainingPtype<-exp(trainingPtype)
trainingPtype[1:4,1:4]
fpkm=read.csv("29.1免疫的表达矩阵.csv",row.names = 1)
group=read.csv("29.2免疫的分组.csv",row.names = 1)

identical(rownames(group),colnames(fpkm))###一定要保证两者是一致的
testExprData = as.matrix(fpkm)
calcPhenotype(trainingExprData = trainingExprData,
              trainingPtype = trainingPtype,
              testExprData = testExprData,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
resultPtype <- read.csv('DrugPredictions.csv', header = T ,stringsAsFactors = F ,check.names = F,row.names = 1)
AA=as.data.frame(colnames(resultPtype))
names(resultPtype)[1] <- "sample"###将第一列重新命名
#LINC01133、FOXD3.AS1、AP001066.1、AP003555.1

resultPtype <- read.csv('DrugPredictions.csv', header = T ,stringsAsFactors = F ,check.names = F)
colnames(resultPtype)[1]="sample"#纵坐标就是每个样本对应药物的评分【sensitive scores】。
###载入自己的分组
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
km <- read.csv("29.2免疫的分组.csv")
km=km[c(1,3)]
colnames(km) <- c("sample","group")
identical(km$sample,resultPtype$sample)###一定要保证一致才能进行合并
rownames(resultPtype)=resultPtype$sample
resultPtype=resultPtype[km$sample,]
library(stringr)
library(dplyr)
fpkm_km_drug <- km %>%
  inner_join(resultPtype)#####inner_join必须要有相同的列名
library(ggplot2)
library(ggpubr)
library(Hmisc)
gsva_mat=fpkm_km_drug
# 批量绘图 --------------------------------------------------------------------
#定义一个基因列表
gene_list <- colnames(gsva_mat)[3:200]####gsva_mat:行名为样本名，列名为要比较的列名

#循环建立箱线图
dir.create("26.7 GDSC_boxplot")#创建一个新的文件夹
for (gene in gene_list) {
  print(gene)
  Figure = ggplot(gsva_mat, aes(x=group, y=gsva_mat[,gene],color=group)) + 
    ylab(gene) + 
    geom_boxplot() +
    geom_point(position=position_jitter(width=0.1)) +
    theme(text = element_text(size=12),#调整字体大小和轴标题的大小
          axis.title = element_text(size=18)) + 
    scale_x_discrete(labels=c("high","low")) + #调整x轴的标签
    theme_minimal(base_size=12, base_family="sans") + scale_color_manual(values = c("#ED545C","#0075C5")) + geom_signif(comparisons=list(c("high","low")), test = "t.test",map_signif_level=TRUE)
  ggsave(filename = paste0(gene, ".pdf"), path = "26.7 GDSC_boxplot",plot = Figure,width = 8,height = 6)
}

#循环建立小提琴图
dir.create("expr_compare_violin")#创建一个新的文件夹
for (gene in gene_list) {
  print(gene)
  Figure = ggplot(expr_use, aes(x=Group, y=expr_use[,gene],fill=Group)) + geom_violin() + geom_boxplot(outlier.colour = "black",width=0.1,fill="white") +
    ylab(gene) +
    theme(text = element_text(size=12),#调整字体大小和轴标题的大小
          axis.title = element_text(size=18)) + 
    scale_x_discrete(labels=c("Normal","Tumor")) + #调整x轴的标签
    theme_minimal(base_size=12, base_family="sans")
  ggsave(filename = paste0(gene, ".pdf"), path = "expr_compare_violin",plot = Figure,width = 8,height = 6)
}

############################################画共同的箱线图
#Acetalax_1804/Afatinib_1032/Axitinib_1021/Buparlisib_1873/Cediranib_1922/Dactinomycin_1911
#Docetaxel_1007/Eg5_9814_1712/Epirubicin_1511/Erlotinib_1168/Foretinib_2040/Gefitinib_1010
#GNE-317_1926/Ipatasertib_1924/Lapatinib_1558/Niraparib_1177/Ribociclib_1632/Sapitinib_1549
#Savolitinib_1936/Temozolomide_1375/Uprosertib_1553

######提取出有意义的药物
fpkm_km_drug=fpkm_km_drug[colnames(fpkm_km_drug)%in%c("Acetalax_1804","Afatinib_1032",
                                                      "Buparlisib_1873","Sapitinib_1549")]

library(ggplot2)
library(tidyr)
library(ggpubr)
fpkm_km_drug$group=km$group
colnames(fpkm_km_drug)=c("Afatinib_sensitivity(IC50)","Sapitinib_sensitivity(IC50)",
                         "Acetalax_sensitivity(IC50)","Buparlisib_sensitivity(IC50)","group")
gene=fpkm_km_drug[,c(5,1:4)]
gene1=gather(gene,variable,value,-group)
#############绘制整体箱线图
pdf("OO1.pdf",width = 22,height = 20)
K=ggplot(gene1, aes(x=group, y=value, fill=group)) +
  geom_boxplot()+
  theme(text = element_text(size=12), # ###调整字体和轴标题大小
        axis.title = element_text(size=18))+
  scale_x_discrete(labels=c("high", "low"))+
  theme_minimal(base_size=12, base_family="sans") + 
  scale_color_manual(values=c("#0075C5", "#ED545C")) +
  geom_signif(comparisons=list(c("high","low")),test="t.test",map_signif_level=TRUE)+
  facet_wrap(~ variable, ncol=21, scales="free_y") ##ncol：按变量划分
##一颗星：P<0.05；二颗星：P<0.01；三颗星：P<0.001

print(K)
dev.off()

library(GSVAdata)
library(GSVA)
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
GOSet1=getGmt("28.10双硫死亡相关基因.gmt")
gene=read.csv("29.1免疫的表达矩阵(取log2).csv",row.names = 1)
range(gene)
GOSetEs=gsva(expr = as.matrix(gene),gset.idx.list =GOSet1,kcdf="Gaussian",parallel.sz=4,method="gsva" )
write.csv(GOSetEs,"28.10双硫死亡分数计算.csv")
getwd()
score=read.csv("13.3全部中位数分组完成与表达矩阵结合.csv")
table(score$group)
score=score[,c(1,1187,1189)]
rownames(score)=score$X
identical(rownames(score),colnames(gene))##必须保证分组是一致的
gene=gene[,rownames(score)]
fe=read.csv("28.10双硫死亡分数计算.csv",row.names = 1)
fe=fe[,rownames(score)]
fe=t(fe)
fe=as.data.frame(fe)
colnames(fe)="score_fe"
identical(rownames(fe),rownames(score))##必须保证分组是一致的
score_all=cbind(fe,score)
score_all=score_all[,-2]
colnames(score_all)=c("ES_disulfidptosis","risk_score","group")
library(ggplot2)
library(ggpubr)
library(ggExtra)
data_use=score_all
colnames(data_use)
target_gene <- 'risk_score'##随机挑选一个基因
target_column <- data_use[,target_gene]
cor_R <- cor(x = target_column,y = data_use[,1],method = 'pearson')#c("pearson", "kendall", "spearman")
cor_P <- cor.test(x = target_column,y = data_use[,1])$p.value
result <- data.frame("target_gene" = character(),
                     "gene_symbol" = character(),
                     "cor_R" = numeric(),
                     "cor_P" = numeric())
gene_list <- colnames(data_use)[1:2]#设置你要进行相关性分析的范围【我们的目的基因可以都放在最后一列，这样子计算的时候就不会包含这一列了】
for (gene in gene_list) {
  print(gene)
  cor_R <- cor(x = target_column,y = data_use[,gene],method = 'pearson')
  cor_P <- cor.test(x = target_column,y = data_use[,gene])$p.value
  temp_result <- data.frame(target_gene = target_gene,
                            gene_symbol = gene,
                            cor_R = cor_R,
                            cor_P = cor_P)
  result <- rbind(result, temp_result)
}
write.csv(result,"28.10双硫死亡相关性结果(无意义).csv")
plot_Example <- ggplot(data_use,aes(x = risk_score,y = ES_disulfidptosis))+
  geom_point(size = 5,color = '#2570A4',alpha = 3)+
  theme_bw()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 12),axis.ticks.length = unit(0.2,'cm'),axis.ticks = element_line(size = 1),panel.border = element_rect(size = 1.5),panel.grid = element_blank())+
  geom_smooth(method = 'lm',se = T,color = 'black',size = 2.0,fill = '#7A8991')+ stat_cor(method = "pearson",digits = 3,size=6)
ggMarginal(plot_Example,type = "density",xparams = list(bw = 0.5, fill = "#B3E283",size = 0.1),yparams = list(bw = 0.5, fill = "#8AB6D6",size = 0.1))


library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(ggplot2)#分面绘图所需
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
# 读取数据：基因差异分析的结果 ----------------------------------------------------------
upGene <- read.csv("1.2DEG_deseq2(mRNAlog2(1.5)).csv",row.names = 1)
table(upGene$regulate)
#upGene=upGene[upGene$regulate!="unchanged",]#GSEA富集分析不需要寻找差异基因
upGene2=upGene[rownames(upGene)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

# GO富集分析与结果存储【GO没有结果】 -------------------------------------------------------------
GO <- enrichGO(gene = Genes$ENTREZID, #输入基因的"ENTREZID"
               OrgDb = org.Hs.eg.db,#注释信息
               keyType = "ENTREZID",
               ont = "all",     #可选条目BP/CC/MF
               pAdjustMethod = "BH", #p值的校正方式
               pvalueCutoff = 0.05,   #pvalue的阈值
               qvalueCutoff = 0.05, #qvalue的阈值
               minGSSize = 5,
               maxGSSize = 5000,
               readable = TRUE)   #是否将entrez id转换为symbol
AA=GO@result
AA=AA[AA$pvalue<0.05,]
write.csv(AA,"29.1 风险分组GO富集分析.csv")

# 结果可视化--柱状图 -----------------------------------------------------------
pdf(file="GO柱状图.pdf",width = 10,height = 8)
barplot(GO, drop = TRUE, 
        showCategory =6,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dev.off()
head(GO@result)

# 结果可视化--点图 ---------------------------------------------------------------
pdf(file="GO点图.pdf",width = 10,height = 8)
dotplot(GO,showCategory = 6,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')  #方式一：分面
dotplot(
  GO,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 15,
  size = NULL,
  split = NULL,
  font.size = 12,
  title = "",
  orderBy = "x",
  label_format = 30)  #方式二：不分面

# KEGG富集分析与结果存储【没有结果】 -----------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)

barplot(KEGG,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 15,
        title = "KEGG_enrichment")#标题
dotplot(KEGG)
View(as.data.frame(KEGG))#以数据框的形式展示KEGG的结果
browseKEGG(KEGG,"hsa04060")#打开KEGG官网查看通路的详细信息

####结果存储
KEGG_results <- DOSE::setReadable(KEGG,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENTREZID") #将ENTREZID转换为genesymbol
AA=KEGG_results@result
AA=AA[AA$pvalue<0.05,]
write.csv(AA,"29.1 风险分组KEGG富集分析.csv")

###############################################################富集分析【重点关注！我们是按照风险的高低分组来进行富集分析的！！】
setwd("D://R/TCGA下载+整理数据/COAD数据/")
exp=read.csv("exp_counts_COAD(全部转化完了).csv",row.names = 1)###最好用全部基因的原因：因为到后面的富集分析就会发现，只有一个基因能够富集到，太过于局限了
##基因过滤
nrow(exp)
#### 常用过滤标准1：仅去除在所有样本里表达量都为零的基因（这是最低的过率标准）
exp1 = exp[rowSums(exp)>0,]##都为0
nrow(exp1)
#### 常用过滤标准2(推荐)：仅保留在一半以上样本里表达的基因——这个条件会导致SLC12A基因的不全
#exp1 = exp1[apply(exp1, 1, function(x) sum(x > 0) > 0.5*ncol(exp1)), ]
###大于或者大于等于，选择大于即可，两者没有明显差异
#nrow(exp1)
############################表达矩阵的处理
table(substr(colnames(exp1),14,16))###查看肿瘤和正常样本的数目
Tumor <- grep('01A',colnames(exp1))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
counts_exp1=exp1[,Tumor]###这里要符合我们的高低风险的分组，所以我们只要肿瘤组织
#表达矩阵以样本为核心去重
library(stringr)
expr_use = counts_exp1[,sort(colnames(counts_exp1))]
k = !duplicated(str_sub(colnames(expr_use),1,12))
table(k)
expr_use = expr_use[,k] 
colnames(expr_use) <- substr(colnames(expr_use),1,12)
range(expr_use)

setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
Group=read.csv("29.2免疫的分组.csv",row.names = 1)
Group1=Group[Group$group=="low",]
Group2=Group[Group$group=="high",]
Group=rbind(Group1,Group2)
count=expr_use[,rownames(Group)]
identical(colnames(count),rownames(Group))#############必须保证一致
write.csv(count,"45.1 用于GSEA富集分析的表达矩阵(全部基因).csv")


####################################开始差异分析
identical(colnames(count),rownames(Group))#############必须保证一致

# 创建分组信息 ------------------------------------------------------------------
table(Group$group_risk.score)
Group <- factor(c(rep("low",times=210),rep("high",times=210)))#创建分组因子变量
Group=relevel(Group,ref = "low")
Data <- data.frame(row.names = colnames(count), #创建分组数据框
                   group = Group)
#基因表达矩阵(counts_exp1)、样本分组数据框(Data)都已经准备完毕


# 开始进行差异表达分析 --------------------------------------------------------------
#DESeq2利用的是count数据
library(DESeq2)
#第一步：构建DEseq2对象(dds)
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = Data,
                              design = ~ group)
#第二步：开始差异分析
dds2 <- DESeq(dds)
res <- results(dds2, contrast=c("group", "high", "low"))#肿瘤在前，对照在后
##或者res= results(dds)
res <- res[order(res$pvalue),]###按照P值排序
summary(res)
my_result <- as.data.frame(res)#转成容易查看的数据框
my_result <- na.omit(my_result)#删除倍数为0的值
#第三步：保存差异分析的结果
library(dplyr)
my_result$Gene_symbol<-rownames(my_result)
my_result <- my_result %>% dplyr::select('Gene_symbol',colnames(my_result)[1:dim(my_result)[2]-1],everything())
rownames(my_result) <- NULL

# DEG的筛选 ------------------------------------------------------------------
my_result$regulate <- ifelse(my_result$padj > 0.05, "unchanged",
                             ifelse(my_result$log2FoldChange > 1, "up-regulated",
                                    ifelse(my_result$log2FoldChange < -1, "down-regulated", "unchanged")))
table(my_result$regulate)
#可以把上调基因和下调基因取出放在一块
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
write.csv(my_result,file= "45.2DEG_deseq2(富集分析1).csv")


########################################################GSEA富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
library(enrichplot)
library(DOSE)
library(pathview)
library(topGO)
#install.packages("ggupset")
library(ggupset)
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
# 制作genelist --------------------------------------------------------------
deg<-read.csv("45.2风险分组_用于富集分析DEG_deseq2(1).csv",header = T,sep = ",")###差异基因
table(deg$regulate)
#deg=deg[deg$regulate!="unchanged",]——GSEA不需要寻找差异基因，只与其原理有关
deg=deg[,-1]
names(deg)[1] <- "SYMBOL"

deg1 <- bitr(geneID = deg$SYMBOL,fromType = "SYMBOL",
             toType = "ENTREZID",OrgDb = org.Hs.eg.db)
DEG <- merge(deg,deg1,by = "SYMBOL")####两个结合是为了有ENTERZID

DEG1<-DEG[,c(1,3,9)]####因为GSEA富集分析是用logFC来进行计算的
colnames(DEG1)
dega<-DEG1[,c(3,2)]###只要后2列
geneList = dega[,2]###提取出logFC一列
names(geneList) = as.character(dega[,1])####进行命名
head(geneList)
geneList = sort(geneList, decreasing = TRUE)###按照从大到小排列
head(geneList)
getwd()
# #########################################gseKEGG分析 ------------------------------------------------------------------
kk2 <- gseKEGG(geneList = geneList,
               organism = 'hsa',
               keyType = "kegg",
               exponent = 1,
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               by = "fgsea")
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
af=as.data.frame(kk2@result)
write.csv(af,"45.3风险分组GSEA富集分析.csv")

# 结果可视化

####画图1：分别取GSEA结果的前5个后5个展示
#num=5
#gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
#gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])

####画图2：取出我想要的
##排序
kk3=kegg_result[(order(kegg_result$enrichmentScore,decreasing = T)),]##decreasing = T:降序排列

pdf(file="GSEA1.pdf",width = 18,height = 16)####如果图片不完整，那么就用宽度来解决
MM=gseaplot2(kk2, geneSetID = rownames(kk3)[c(1,2,3,4,7,49,48,47,46,45)])
print(MM)
dev.off()


pdf(file="GSEA3.pdf",width = 18,height = 16)####如果图片不完整，那么就用宽度来解决
MM=gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
print(MM)
dev.off()


# 单独展示某一个条目 ---------------------------------------------------------------
pdf(file="GSEA4.pdf",width = 18,height = 16)####如果图片不完整，那么就用宽度来解决
P=gseaplot2(kk2,
            title = "Ferroptosis",  #设置标题
            "hsa04216", #绘制hsa04658通路的结果，通路名称与编号对应
            color="red", #线条颜色
            base_size = 20, #基础字体的大小
            subplots = 1:3, 
            pvalue_table = T) # 显示p值
print(P)
dev.off()
# GSEA 结果表达分布的脊线图 ---------------------------------------------------------
#跑分和预排序是可视化 GSEA 结果的传统方法,可视化基因集的分布和富集分数
ridgeplot(kk2,
          showCategory = 50,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 32)
gseaplot(kk2, geneSetID = 1, by = "all", title = kk2$Description[1])




#####################################gsego分析 -----------------------------------------------------------------
gsea_go <- gseGO(geneList     = geneList,#根据LogFC排序后的基因列表
                 OrgDb        = org.Hs.eg.db,
                 ont          = "ALL",#GO分析的模块
                 minGSSize    = 10,#最小基因集的基因数
                 maxGSSize    = 500,#最大基因集的基因数
                 pvalueCutoff = 0.05,#p值的阈值
                 verbose      = FALSE)#是否输出提示信息

class(gsea_go)
colnames(gsea_go@result)
kegg_result <- as.data.frame(gsea_go)
rownames(gsea_go@result)[head(order(gsea_go@result$enrichmentScore))]
af=as.data.frame(gsea_go@result)
write.csv(af,"GSEA——GO富集分析.csv")
###先看排序
kk3=af[(order(af$enrichmentScore)),]
View(as.data.frame(gsea_go))
#绘图
gseaplot2(gsea_go, geneSetID = rownames(gsea_go@result)[head(order(gsea_go@result$enrichmentScore),num)])


# 单独展示某一个条目 ---------------------------------------------------------------
gseaplot2(gsea_go,
          title = "regulation of lymphocyte activation",  #设置标题
          "GO:0051249", #绘制hsa04658通路的结果，通路名称与编号对应
          color="red", #线条颜色
          base_size = 20, #基础字体的大小
          subplots = 1:3, 
          pvalue_table = T) # 显示p值

##################################################################利用MSigDB是GSEA官方提供的基因集数据库，整合了数万个注释基因集资源
#devtools::install_github("ToledoEM/msigdf")
library(msigdf)
library(dplyr)
##以C7为例【C7是免疫相关的基因集合】：
#提取C7注释(human)：
c7 <- msigdf.human %>%
  filter(category_code == "c7") %>% select(geneset, symbol) %>% as.data.frame
head(c7)
#准备genelist文件：
#这里genelist需要的是symbol号；
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
deg<-read.csv("45.2风险分组_用于富集分析DEG_deseq2(1).csv",row.names = 1)###差异基因
table(deg$regulate)
genelist2 <- deg$log2FoldChange
names(genelist2) <- deg$Gene_symbol
#添加entrez ID列：
##symbol转entrez ID：
symbol <- deg$Gene_symbol
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)
#genelist缺失基因过滤(ID转换中缺失的部分)
genelist2 <- genelist2[names(genelist2) %in% entrez[,1]]
length(genelist2)
head(genelist2)
#将genelist按照log2FC值从高到低进行排序：
genelist2 <- sort(genelist2, decreasing = T)
head(genelist2)
c2_ges <- GSEA(genelist2,
               TERM2GENE = c7,
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose = FALSE,
               eps = 0)
c2_ges_result <- c2_ges@result
write.csv(c2_ges_result, file = c('45.6 免疫——GSEA(MSigDB-c7).csv'))


# 结果可视化--分别取GSEA结果的前5个后5个展示 -----------------------------------------------
#View(as.data.frame(kk2@result))
#先排个序看一下
kk3=c2_ges_result[(order(c2_ges_result$enrichmentScore,decreasing = T)),]##decreasing = T:降序排列

pdf(file="GSEA6.pdf",width = 18,height = 16)####如果图片不完整，那么就用宽度来解决
#按照以下的方式模仿
#num=5
#gseaplot2(gsea_go, geneSetID = rownames(c2_ges@result)[head(order(c2_ges@result$enrichmentScore),num)])

MM=gseaplot2(c2_ges, geneSetID = rownames(kk3)[c(10,12,13,17,20)])
print(MM)
dev.off()


# 单独展示某一个条目 ---------------------------------------------------------------
pdf(file="GSEA4.pdf",width = 18,height = 16)####如果图片不完整，那么就用宽度来解决
P=gseaplot2(kk2,
            title = "Ferroptosis",  #设置标题
            "hsa04216", #绘制hsa04658通路的结果，通路名称与编号对应
            color="red", #线条颜色
            base_size = 20, #基础字体的大小
            subplots = 1:3, 
            pvalue_table = T) # 显示p值
print(P)
dev.off()
# GSEA 结果表达分布的脊线图 ---------------------------------------------------------
#跑分和预排序是可视化 GSEA 结果的传统方法,可视化基因集的分布和富集分数
ridgeplot(kk2,
          showCategory = 50,
          fill = "p.adjust",
          core_enrichment = TRUE,
          label_format = 32)
gseaplot(kk2, geneSetID = 1, by = "all", title = kk2$Description[1])





####################################################################cibersort免疫浸润分析
setwd("D://R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
library(CIBERSORT)
library(ggplot2)#绘
library(pheatmap)#绘制热图
library(ggpubr)#堆积比例图
library(reshape2)#数据处理
library(tidyverse)#数据处理
# 免疫细胞类型列表 ----------------------------------------------------------------
fpkm=read.csv("29.3免疫的表达矩阵(未取log2).csv",row.names = 1)###不需要取log2
LM22_local <- read.table("LM22.txt",header = T,row.names = 1,sep = "\t")
data(LM22)#调用cibersort包中自带的LM22数据，22种免疫细胞参考marker基因表达情况
all(LM22==LM22_local)#判断官网下载的文件是否与cibersort包中自带的LM22数据一致
# 正式开始Cibersort免疫细胞浸润分析 ---------------------------------------------------
#运行CIBERSORT以估计样本中免疫细胞类型的比例
result <- cibersort(sig_matrix = LM22, mixture_file =fpkm , perm = 1000, QN = F)
#参数包括 sig_matrix、mixture_file、perm 和 QN。
#sig_matrix 参数是 CIBERSORT 软件包的内置数据，mixture_file 参数是待测样本的基因表达数据。perm 参数表示是否使用随机排列法，QN 参数是 TRUE 或 FALSE，表示是否使用质量归一化(芯片数据设置为T，测序数据就设置为F）。cibersort 函数会返回一个数据框，包含了各种免疫细胞类型的比例和基因表达数据的质量归一化结果。
#一般perm设置是1000
result <- as.data.frame(result)
colnames(result)
write.csv(result,"32.1cibersort_result(1000).csv")

#`P-value`越小越可信；
#Correlation原表达矩阵乘以细胞占比后的数据矩阵与原表达矩阵的相关性
#RMSE 均方根误差，越小效果越好

# 结果可视化 -------------------------------------------------------------------
getwd()
result=read.csv("32.1cibersort_result(1000).csv",row.names = 1)
result1 <- result[,1:ncol(LM22)]##指1：22种免疫细胞
result1 <- result1[,apply(result1, 2, function(x){sum(x)>0})]#删除全是0的列
#在矩阵的列上应用函数：apply(X, 2, fun)，其中X是矩阵，fun是要应用的函数，2表示按列应用函数。
#热图
pheatmap(result1,
         color = colorRampPalette(c("#4CB43C", "#FEFCFB", "#ED5467"))(100),
         border="skyblue",#边框颜色
         main = "Heatmap",#指定图表的标题
         show_rownames = F,#是否展示行名
         show_colnames = T,#是否展示列名
         cexCol = 1,#指定列标签的缩放比例。
         scale = 'row',#指定是否应按行方向或列方向居中和缩放，或不居中和缩放。对应的值为row, column和none。
         cluster_col=T,#分别指定是否按列和行聚类。
         cluster_row=F,
         angle_col = "45",#指定列标签的角度。
         legend = F,#指定是否显示图例。
         legend_breaks=c(-3,0,3),#指定图例中显示的数据范围为-3到3。
         fontsize_row = 10,#分别指定行标签和列标签的字体大小。
         fontsize_col = 10)
dev.off()


# 盛夏的果实函数科普一： --------------------------------------------------------------
####在R语言中，melt函数是用来将一个数据集转化为“长格式”的。##就是把宽转换为长
###所谓的“长格式”是指将一个数据集的多个变量的值放在一个列中，
###而不是将每个变量放在不同的列中。这样就可以使用一列来记录变量的名称，
###另一列来记录变量的值。

# ggplot2绘制堆积比例图 ----------------------------------------------------------
#画堆积图的原因：每个样本的22种免疫细胞的比例相加为1
#数据整理
Data=read.csv("29.2免疫的分组.csv",row.names = 1)
result1=result1[rownames(Data),]
identical(rownames(result1),rownames(Data))
data <- cbind(rownames(result1),result1)
colnames(data)[1] <- "Samples"
data <- melt(data,id.vars = c("Samples"))###将所有免疫细胞类型作为一列
colnames(data) <- c('Samples','celltype','proportion')
#开始绘图
mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C',"#ACE89C")
ggplot(data,
       aes(Samples,proportion,fill=celltype))+geom_bar(stat="identity",position="fill")+#x 轴是变量 Samples，y 轴是变量 proportion，条形的填充颜色由变量 celltype 决定
  scale_fill_manual(values=mycolors)+#填入需要填充的颜色，22种免疫细胞
  ggtitle("Proportion of immune cells")+theme_gray()+theme(axis.ticks.length=unit(3,'mm'),axis.title.x=element_text(size=11))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
  guides(fill=guide_legend(title="Types of immune cells"))



# 盛夏的果实函数科普二 --------------------------------------------------------------
#pivot_longer是tidyr包中的一个函数，它可以将一个数据框中的多个列转换为两列：
#一列是原来的列名，另一列是对应的值。
#pivot_longer(data, cols, names_to = "新的列名", values_to = "新的值列名")

# ggplot2绘制分组箱式图 ----------------------------------------------------------
#【这个图经常在文献中看到的】

# 提取其中一组（肿瘤/正常）绘图 ---------------------------------------------------------
data2 <- cbind(result1,Data)
data2 <- data2[data2$group_risk.score=='low',]
data2 <- pivot_longer(data = data2,
                      cols = 1:22,
                      names_to = "celltype",
                      values_to = "proportion")
sum(data2$proportion)
#pdf(file="单组箱式图.pdf",width = 10,height = 8)
ggboxplot(data = data2,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",
          color = "black",
          xlab = "Types of immune cells",#x 轴标签。
          ylab = NULL,
          title = "TME Cell composition",
          fill = "celltype",
          legend.position = "bottom",
          ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) 
#dev.off() 


# ggplot2绘制分组箱式图升级版+显著性P值 -------------------------------------------------
data3 <- cbind(result1,Data)
colnames(data3)
data3 <- data3[,c(23,24,1:22)]
data3 <- pivot_longer(data = data3,
                      cols = 3:24,
                      names_to = "celltype",
                      values_to = "proportion")
#开始绘图
#pdf(file="分组箱式图+显著性P值.pdf",width = 10,height = 8)
ggboxplot(data = data3,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",#绘制箱形图的响应变量。
          combine = TRUE,#是否将数据合并为一个箱形图。
          merge = FALSE,#是否将相同值的分组合并。
          color = "black",#箱形图边框的颜色。
          fill = "group",#箱形图填充色。
          palette = c("#ED5462","#81D5B0"),#颜色调色板。
          title = "TME Cell composition",#图形标题。
          xlab = NULL,#x 轴标签。
          ylab = "Cell composition",#y 轴标签
          bxp.errorbar = FALSE,#是否在箱形图中绘制误差条。
          bxp.errorbar.width = 0.2,#误差条宽度。
          facet.by = NULL,#基于哪些变量进行分面
          panel.labs = NULL,#分面的标签
          short.panel.labs = TRUE,#是否将分面标签缩短
          linetype = "solid",#线条类型
          size = NULL,#图形大小。
          width = 0.8,#箱形图的宽度。
          notch = FALSE,#是否在箱形图中绘制刻度。
          outlier.shape = 20,#异常值标记的形状。
          select = NULL,#要绘制的变量
          remove = NULL,#不要绘制的变量。
          order = NULL,#箱形图的排序方式。
          error.plot = "pointrange",#如何绘制误差，可以是 "pointrange" 或 "errorbar"。
          label = NULL,#要添加的标签
          font.label = list(size = 12, color = "black"),#标签的字体属性
          label.select = NULL,#要添加标签的数据点
          repel = TRUE,#是否使用 repel 库的功能使标签互不重叠
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) +  #这个函数的作用是将 x 轴文本旋转 90 度，并调整其对齐方式。
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.",hide.ns = T) 
#dev.off()  
data4=data3[data3$celltype=="Dendritic.cells.activated",]
data5=data4[data4$group=="low",]
mean(data5$proportion)

ggboxplot(data = data4,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",
          combine = TRUE,
          merge = FALSE,#是否将相同值的分组合并。
          color = "black",#箱形图边框的颜色。
          fill = "group",#箱形图填充色。
          palette = c("#ED5462","#81D5B0"))

###############################################免疫细胞的相关性热图
# 加载需要的R包 -----------------------------------------------------------------
library(ggcorrplot)
library(tidyr)

# 导入数据 --------------------------------------------------------------------
DEG_expr <- read.csv("6.临床数据与PANoptosis—lncRNAs表达数据结合.csv",row.names = 1)#行名基因名，列名样本名
DEG_expr=DEG_expr[,-c(1:9)]
DEG_expr=DEG_expr[rownames(Data),]
range(DEG_expr)
data=DEG_expr
# cibersort运算结果【用的是这个结果来进行演示的】
cibersort_result <- read.csv("32.1cibersort_result(1000).csv",row.names = 1)#行名样本名，列名免疫细胞浸润评分
#数据格式整理
#data <- t(DEG_expr)#数据转置，行名是样本名，列名基因名
cibersort_result <- cibersort_result[,1:22]##提取出免疫细胞亚型

# 解决历史遗留问题 ----------------------------------------------------------------
#问题一：对于Cibersort免疫细胞浸润分析，如果输入的基因表达矩阵是以log转换的形式表示的，
#需要进行反转换以获得原始表达值，因为Cibersort需要原始的表达值来进行细胞比例估计。

#进行log逆转换(此处以2为底举例，实际以GEO官网数据集介绍为准)
#if(max(data) < 50) {data <- 2^data}【这个就不需要了】

#问题二：在Cibersort中，输入基因表达矩阵应该是进行标准化处理后的数据。
#如果对数转换后再进行反转换，得到的数据是未标准化的原始数据，这与Cibersort需要的标准化数据是不一样的。
#data <- t(DEG_expr)#数据转置，行名是样本名，列名基因名
#check一下
identical(rownames(data),rownames(cibersort_result))

#数据简单处理
# Shapiro-Wilk检验#【相关性分析——检验正态分布数据】
shapiro.test(data$SHISA3)#适用于小样本数据，N≤50
# Kolmogorov-Smirnov检验(柯尔莫哥洛夫检验)
ks.test(data$LINC01614, pnorm, mean = mean(data$LINC01614), sd = sd(data$LINC01614))#适用于大样本数据，N＞50

# 计算相关性 -------------------------------------------------------------------
colnames(data)
colnames(cibersort_result)
target_Gene <- "AL365181.3"#确定目标基因
expr <- as.numeric(data[,target_Gene])
#举例：计算SHISA3和T.cells.CD8的相关性分析
cor.test(cibersort_result$T.cells.CD8, expr, use = "everything", method = "spearman")

# 封装一个函数来快速计算 -------------------------------------------------------------
# spearman相关性分析
calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "spearman")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}
# pearson相关性分析【我们这里用的是pearson相关性分析】
calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "pearson")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}

# 计算4个基因与免疫细胞浸润分数的相关性 -----------------------------------------------------
# 选择要计算相关性的列[我们是计算相关的列]
data_to_calculate <- data[,colnames(data)%in%c("LINC01133","FOXD3.AS1","AP001066.1","AP003555.1")]
data_to_calculate=as.data.frame(data_to_calculate)

# 新建一个空的数据框保存结果
results <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
# 使用for循环遍历数据框中的每一列，并计算相关性
for (i in 1:ncol(data_to_calculate)) {
  print(i)
  gene_expr <- data_to_calculate[, i]
  corr_result <- calculate_correlation(gene_expr, cibersort_result)###这个函数是通过上面的自定义函数，自己来选择
  # 将每次计算的结果添加到新的数据框中
  results <- rbind(results, corr_result)
}
# 查看列名，下方修改的依据
colnames(data_to_calculate)

# 手动修改每一个基因名称
results$Gene <- c(rep("LINC01133", 22),
                  rep("FOXD3.AS1", 22),
                  rep("AP001066.1", 22),
                  rep("AP003555.1", 22))
write.csv(results,"32.5计算4个基因与免疫细胞浸润分数的相关性.csv")
results=results[results$p.value<0.05,]##一定要挑出来这P<0.05，才能够后续画图
results <- results[,c(1,2,3)]#转换数据，这里不需要显著性p值，也可以通过显著性P值进行筛选
colnames(results)
# 使用spread函数将长数据转换为宽格式
results_wide <- spread(results, key = im_cell, value = Cor)
# 将基因名设置为行名
rownames(results_wide) <- results_wide$Gene
# 删除原始数据框中的基因名列
results_wide$Gene <- NULL

# 绘制相关性矩阵图
ggcorrplot(t(results_wide), 
           # 显示颜色图例
           show.legend = T, 
           # 设置相关系数矩阵的颜色，其中第一个为最小值对应的颜色，
           # 中间的白色为0对应的颜色，最后一个为最大值对应的颜色
           colors = c("#2166AC", "white", "#B2182B"), 
           # 设置数字显示的位数
           digits = 2, 
           # 显示变量标签（默认为TRUE）
           lab = T)





setwd("D://R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
gene=read.csv("29.1免疫的表达矩阵(取log2).csv",row.names = 1)
gene=gene[rownames(gene)%in%c("BTLA","HHLA2","IDO2"),]####一定要基于前人的文献才可以进行选择免疫检查点
################行名为样本，列名为基因名
gene=t(gene)
gene=as.data.frame(gene)
group=read.csv("29.2免疫的分组.csv",row.names = 1)
identical(rownames(gene),rownames(group))
gene=gene[rownames(group),]
gene$group=group$group
gene=gene[c(4,1:3)]#####保证group在最前面
write.csv(gene,"42.2免疫检查点相关性矩阵.csv")

gene=read.csv("42.2免疫检查点相关性矩阵.csv",row.names = 1)
#################################画图
library(ggplot2)
library(tidyr)
library(ggpubr)

############################################将数据框从宽格式转变为长格式
gene1=gather(gene,variable,value,-group)
#############绘制整体箱线图
pdf("OO33.pdf",width = 22,height = 20)
K=ggplot(gene1, aes(x=group, y=value, fill=group)) +
  geom_boxplot()+
  theme(text = element_text(size=12), # ###调整字体和轴标题大小
        axis.title = element_text(size=18))+
  scale_x_discrete(labels=c("high", "low"))+
  theme_minimal(base_size=12, base_family="sans") + 
  scale_color_manual(values=c("#0075C5", "#ED545C")) +
  geom_signif(comparisons=list(c("high","low")),test="t.test",map_signif_level=TRUE)+
  facet_wrap(~ variable, ncol=21, scales="free_y") ##ncol：按变量划分
##一颗星：P<0.05；二颗星：P<0.01；三颗星：P<0.001

print(K)
dev.off()

##################################################3相关性分析


#################################批量计算基因间的相关性
# 计算相关性R和显著性p值 ------------------------------------------------------------
data_use=gene
target_gene <- c("score")
#################################批量计算基因间的相关性
#生成一个空白的数据框，存放结果
result <- c()
#创建要进行相关性分析的基因列表
gene_list <- colnames(data_use)[2:22]#进行相关性分析的范围【我们是选择全部的基因】
#循环读取
for (gene in gene_list) {
  for(target in target_gene){
    message(gene,target)#打印【报出信息】
    cor_R <- cor.test(x = data_use[,target],y = data_use[,gene],method = 'pearson')###基因近似正态分布
    P.value<-cor_R[["p.value"]]
    cor=cor_R[["estimate"]][["cor"]]
    temp_result<-c(gene,target,cor,P.value)#让这个变量一直在变
    result <- rbind(result, temp_result)
  }
  
}
colnames(result)=c("mRNA","score","cor_R","cor_P")
write.csv(result,"42.2 免疫检查点相关性的分析.csv")
################################################提取出4个预后基因的的fpkm_lncRNA矩阵与临床数据的结合
setwd("D://R/TCGA下载+整理数据/COAD数据/")
exp=read.csv("output_mRNA_lncRNA_expr/TCGA-COAD_lncrna_expr_fpkm.csv",row.names = 1)
############################fpkm表达矩阵的处理
table(substr(colnames(exp),14,16))###查看肿瘤和正常样本的数目
Tumor <- grep('01A',colnames(exp))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_exp=exp[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
Normal <- grep('11A',colnames(exp))
Normal #正常样本所处位置
Normal_exp=exp[,Normal]
fpkm_exp=cbind(Normal_exp,Tumor_exp)
######提取出fpkm矩阵中的相关lncRNA
fpkm_exp1=fpkm_exp[rownames(fpkm_exp)%in%c("LINC01133","FOXD3-AS1","AP001066.1","AP003555.1"),]
############不要忘了log2的放缩
range(fpkm_exp1)
fpkm_exp1=log2(fpkm_exp1+1)###这里已经放缩了


###################临床信息与表达矩阵的结合
#表达矩阵
table(substr(colnames(fpkm_exp1),14,16))#查看肿瘤样本和正常样本数量
##只保留肿瘤样本，把正常的样本剔除掉
Tumor <- grep('01A',colnames(fpkm_exp1))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_lnc=fpkm_exp1[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
library(stringr)
expr_use = Tumor_lnc[,sort(colnames(Tumor_lnc))]
k = !duplicated(str_sub(colnames(expr_use),1,12))
table(k)
expr_use = expr_use[,k] 
colnames(expr_use) <- substr(colnames(expr_use),1,12)
##临床信息
setwd("D://R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
survival_data=read.csv("1.3结肠癌临床数据.csv",row.names = 1)
survival_data=survival_data[c(1,5,2,3,4,6,7,8,9)]
table(substr(rownames(survival_data),14,16))#查看肿瘤样本和正常样本数量
######只要肿瘤数据
Tumor <- grep('01A',rownames(survival_data))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
survival_data=survival_data[Tumor,]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
survival_data = survival_data[sort(rownames(survival_data)),]
library(stringr)
k = !duplicated(str_sub(rownames(survival_data),1,12))
table(k)
survival_data = survival_data[k,] 
rownames(survival_data) <- substr(rownames(survival_data),1,12)
library(stringr)
rownames(survival_data)=str_replace_all(rownames(survival_data),"-",".")
#将表达矩阵和生存信息合并
datExpr <- t(expr_use)
datExpr=as.data.frame(datExpr)
dim(datExpr)
dim(survival_data)
head(rownames(datExpr))
head(rownames(survival_data))
s = intersect(rownames(datExpr),rownames(survival_data));length(s)
####通过取交集，进行匹配【有的数据是只有表达数据，却没有临床信息，这样的样本是没有意义的】
datExpr1 = datExpr[s,]
survival_data2 = survival_data[s,]
dim(datExpr1)
dim(survival_data2)
identical(rownames(survival_data2),rownames(datExpr1))
#生存信息与表达矩阵合并
datExpr <- cbind(survival_data2,datExpr1)####只有这样才能够保持time和event在最前面
#保存数据
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
write.csv(datExpr,"44.1 临床数据与4个预后基因lncRNAs表达数据结合.csv")


################################################mRNA_FPKM表达矩阵与临床数据的结合
setwd("D://R/TCGA下载+整理数据/COAD数据/")
exp=read.csv("output_mRNA_lncRNA_expr/TCGA-COAD_mrna_expr_fpkm.csv",row.names = 1)
#### 过滤标准1：仅去除在所有样本里表达量都为零的基因（这是最低的过率标准）
exp = exp[rowSums(exp)>0,]##都为0
nrow(exp)
############################fpkm表达矩阵的处理
table(substr(colnames(exp),14,16))###查看肿瘤和正常样本的数目
Tumor <- grep('01A',colnames(exp))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_exp=exp[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
Normal <- grep('11A',colnames(exp))
Normal #正常样本所处位置
Normal_exp=exp[,Normal]
fpkm_exp1=cbind(Normal_exp,Tumor_exp)
############不要忘了log2的放缩
range(fpkm_exp1)
fpkm_exp1=log2(fpkm_exp1+1)###这里已经放缩了


###################临床信息与表达矩阵的结合
#表达矩阵
table(substr(colnames(fpkm_exp1),14,16))#查看肿瘤样本和正常样本数量
##只保留肿瘤样本，把正常的样本剔除掉
Tumor <- grep('01A',colnames(fpkm_exp1))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_lnc=fpkm_exp1[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
library(stringr)
expr_use = Tumor_lnc[,sort(colnames(Tumor_lnc))]
k = !duplicated(str_sub(colnames(expr_use),1,12))
table(k)
expr_use = expr_use[,k] 
colnames(expr_use) <- substr(colnames(expr_use),1,12)
##临床信息
setwd("D://R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
survival_data=read.csv("1.3结肠癌临床数据.csv",row.names = 1)
survival_data=survival_data[c(1,5,2,3,4,6,7,8,9)]
table(substr(rownames(survival_data),14,16))#查看肿瘤样本和正常样本数量
######只要肿瘤数据
Tumor <- grep('01A',rownames(survival_data))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
survival_data=survival_data[Tumor,]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
survival_data = survival_data[sort(rownames(survival_data)),]
library(stringr)
k = !duplicated(str_sub(rownames(survival_data),1,12))
table(k)
survival_data = survival_data[k,] 
rownames(survival_data) <- substr(rownames(survival_data),1,12)
library(stringr)
rownames(survival_data)=str_replace_all(rownames(survival_data),"-",".")
#将表达矩阵和生存信息合并
datExpr <- t(expr_use)
datExpr=as.data.frame(datExpr)
dim(datExpr)
dim(survival_data)
head(rownames(datExpr))
head(rownames(survival_data))
s = intersect(rownames(datExpr),rownames(survival_data));length(s)
####通过取交集，进行匹配【有的数据是只有表达数据，却没有临床信息，这样的样本是没有意义的】
datExpr1 = datExpr[s,]
survival_data2 = survival_data[s,]
dim(datExpr1)
dim(survival_data2)
identical(rownames(survival_data2),rownames(datExpr1))
#生存信息与表达矩阵合并
datExpr <- cbind(survival_data2,datExpr1)####只有这样才能够保持time和event在最前面
#保存数据
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
write.csv(datExpr,"44.2 临床数据与mRNAs表达数据结合.csv")


############################################################lncRNA与mRNA的相关性分析
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
exp1=read.csv("44.2 临床数据与mRNAs表达数据结合.csv",row.names = 1)
exp1=exp1[,-c(1:9)]
exp2=read.csv("44.1 临床数据与4个预后基因lncRNAs表达数据结合.csv",row.names = 1)
exp2=exp2[,-c(1:9)]

########################################矩阵的合并
exp5=cbind(exp2,exp1)
range(exp5)####计算相关性可以暂时不用log2缩小差距

#################################批量计算基因间的相关性
# 计算相关性R和显著性p值 ------------------------------------------------------------
data_use=exp5
#"LINC01133","FOXD3-AS1","AP001066.1","AP003555.1"
target_gene <- c("LINC01133")

#################################批量计算基因间的相关性
#生成一个空白的数据框，存放结果
result <- c()
#创建要进行相关性分析的基因列表
gene_list <- colnames(data_use)[1:19490]#进行相关性分析的范围【我们是选择全部的基因】
#循环读取
for (gene in gene_list) {
  for(target in target_gene){
    message(gene,target)#打印【报出信息】
    cor_R <- cor.test(x = data_use[,target],y = data_use[,gene],method = 'pearson')###基因近似正态分布
    P.value<-cor_R[["p.value"]]
    cor=cor_R[["estimate"]][["cor"]]
    temp_result<-c(gene,target,cor,P.value)#让这个变量一直在变
    result <- rbind(result, temp_result)
  }
  
}
colnames(result)=c("mRNA","LINC01133","cor_R","cor_P")
result=as.data.frame(result)
result=na.omit(result)###会出现缺失值
result$cor_R=as.numeric(result$cor_R)
result$cor_P=as.numeric(result$cor_P)
#自行设置阈值（二选一），筛选要的基因对
table(result$cor_R >0.3)#查看相关性R大于0.3的基因个数
table(result$cor_P < 0.05)#查看p值小于0.05（阈值自选）的基因个数
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
write.csv(result,"44.3LINC01133_microRNA的分析.csv")

index <- result$cor_P < 0.05
index1<-result$cor_R >0.3#注意：正常情况——相关性系数在【-1，1】之间都是相关性比较高的，但是我们这里选择的是与基因正相关的
result_filt <- result[index&index1,]
########去重复
library(dplyr)
result_filt1 <- result_filt[!duplicated(result_filt$mRNA),]
rownames(result_filt1)=result_filt1$mcRNA
write.csv(result_filt1,"44.3LINC01133_microRNA的分析(0.1).csv")


################################################################4个蛋白取交集
L1=read.csv("44.3LINC01133_mRNA的分析(0.3).csv",row.names = 1)
colnames(L1)=c("mRNA","lncRNA","cor_R","cor_P")
L2=read.csv("44.4FOXD3-AS1_mRNA的分析(0.3).csv",row.names = 1)
colnames(L2)=c("mRNA","lncRNA","cor_R","cor_P")
L3=read.csv("44.5AP001066.1_mRNA的分析(0.3).csv",row.names = 1)
colnames(L3)=c("mRNA","lncRNA","cor_R","cor_P")
L4=read.csv("44.6AP003555.1_mRNA的分析(0.3).csv",row.names = 1)
colnames(L4)=c("mRNA","lncRNA","cor_R","cor_P")
L=rbind(L1,L2,L3,L4)
L=L[!duplicated(L$mRNA),]
write.csv(L,"44.7LINC01133_mRNA的相关性（去重复）.csv")




############################富集分析 ------------------------------------------------------------
library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(ggplot2)#分面绘图所需
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
# 读取数据：基因差异分析的结果 ----------------------------------------------------------
upGene1 <- read.csv("54.1 hsa.mir.17_mRNA的分析(交集).csv",row.names = 1)
upGene2 <- read.csv("54.2 hsa.mir.93_mRNA的分析(交集).csv",row.names = 1)
colnames(upGene1)=c("mRNA","miRNA","cor_R","cor_P")
colnames(upGene2)=c("mRNA","miRNA","cor_R","cor_P")
upGene=rbind(upGene1,upGene2)
upGene=upGene[!duplicated(upGene$mRNA),]
rownames(upGene)=upGene$mRNA
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

# GO富集分析与结果存储【GO没有结果】 -------------------------------------------------------------
GO <- enrichGO(gene = Genes$ENTREZID, #输入基因的"ENTREZID"
               OrgDb = org.Hs.eg.db,#注释信息
               keyType = "ENTREZID",
               ont = "all",     #可选条目BP/CC/MF
               pAdjustMethod = "BH", #p值的校正方式
               pvalueCutoff = 0.05,   #pvalue的阈值
               qvalueCutoff = 0.05, #qvalue的阈值
               minGSSize = 5,
               maxGSSize = 5000,
               readable = TRUE)   #是否将entrez id转换为symbol

write.csv(GO_result,"55.5 lncRNA_microRNA_mRNA的GO富集分析.csv")

# 基础柱状图 -------------------------------------------------------------------
#barplot(GO)
#barplot(GO, drop = TRUE, 
#        showCategory =6,split="ONTOLOGY") + 
#  facet_grid(ONTOLOGY~., scale='free')


# 或许你想展示某些特定条目 ------------------------------------------------------------
GO_result <- as.data.frame(GO)
GO_result=GO_result[c(1,4,68,69,70),]
GO_result$logPvalue <- (-log(GO_result$pvalue))
GO_result=arrange(GO_result,logPvalue)
ggplot(GO_result, aes(x = Description, y = logPvalue, fill = Description)) + #指定数据为merge_data，以Description作为x轴，以logPvalue作为y轴，以ONTOLOGY作为颜色填充。
  geom_bar(stat = "identity", width = 0.8, color = "black") + #添加一个柱状图层，使用identity统计方法（不进行任何变换），每个柱子的宽度为0.8，颜色为黑色
  coord_flip() + #翻转坐标轴，交换x轴和y轴，使得柱状图变为横向。
  scale_y_continuous(expand = c(0,0))+ #设置y轴的范围为连续型变量（continuous），设置expand参数为c(0,0)，意味着y轴不需要扩展（不需要增加额外的空白）。
  theme(panel.background = element_rect(fill = "white"))




####################单独想提取出我想要的那几行
BP_top10<- GO_result[c(15,80,96),]   # 筛选 ONTOLOGY 为 BP 的子集
BP_top10%>%arrange(BP_top10$pvalue) 


MF_top10 <- GO_result %>%
  filter(ONTOLOGY == "MF") %>%  # 筛选 ONTOLOGY 为 MF 的子集
  arrange(pvalue) %>%        # 按照 p.adjust 的升序排列
  head(3) 

CC_top10 <- GO_result %>%
  filter(ONTOLOGY == "CC") %>%  # 筛选 ONTOLOGY 为 MF 的子集
  arrange(pvalue) %>%        # 按照 p.adjust 的升序排列
  head(3)                      # 提取前5行

merge_data <- rbind(BP_top10,CC_top10,MF_top10)
merge_data$ONTOLOGY <- factor(merge_data$ONTOLOGY, levels = c("BP", "CC","MF"))
merge_data$logPvalue <- (-log(merge_data$pvalue))
# 美化版本 --------------------------------------------------------------------
colnames(merge_data)
merge_data$Description <- factor(merge_data$Description, levels = unique(merge_data$Description[order(merge_data$ONTOLOGY)]))

ggplot(merge_data, aes(x = Description, y = logPvalue, fill = ONTOLOGY)) + #指定数据为merge_data，以Description作为x轴，以logPvalue作为y轴，以ONTOLOGY作为颜色填充。
  geom_bar(stat = "identity", width = 0.8, color = "black") + #添加一个柱状图层，使用identity统计方法（不进行任何变换），每个柱子的宽度为0.8，颜色为黑色
  scale_x_discrete(limits = unique(merge_data$Description[order(merge_data$ONTOLOGY)])) + #设置x轴的标签顺序，根据ONTOLOGY对Description排序，以便按照ONTOLOGY的顺序显示。unique()函数保证每个标签只出现一次。
  coord_flip() + #翻转坐标轴，交换x轴和y轴，使得柱状图变为横向。
  scale_y_continuous(expand = c(0,0))+ #设置y轴的范围为连续型变量（continuous），设置expand参数为c(0,0)，意味着y轴不需要扩展（不需要增加额外的空白）。
  theme(panel.background = element_rect(fill = "white"))+ #设置背景为白色，使用element_rect()函数设置元素的矩形形状
  theme(axis.text = element_text(size = 16, family = "Arial"),
        axis.title = element_text(size=18,family="Arial",colour = "black"),
        legend.text = element_text(size = 16, family = "Arial"),
        legend.title = element_text(size = 18, family = "Arial"),
        legend.position = "top")+ #设置各种文本元素的字体、大小和颜色，包括坐标轴文本、坐标轴标题、图例文本和图例标题。将图例放置在顶部。
  labs(x = "GO term", y = "-log P-value")+
  scale_fill_nejm(alpha = 0.8) #使用NEJM（New England Journal of Medicine）颜色调色板，设置颜色透明度为0.8，对填充颜色进行比例缩放。


# KEGG富集分析与结果存储【没有结果】 -----------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

####结果存储
KEGG_results <- DOSE::setReadable(KEGG,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENTREZID") #将ENTREZID转换为genesymbol
aa=KEGG@result
aa=aa[aa$pvalue<0.05,]
#KEGG_results <- as.data.frame(KEGG)
write.csv(KEGG_results,"55.4 lncRNA_microRNA_mRNA的KEGG富集分析.csv")

#################################################基础绘图
barplot(KEGG,
        x = "GeneRatio",
        color = "pvalue",
        showCategory = 15,
        title = "KEGG_enrichment")#标题
dotplot(KEGG)
View(as.data.frame(KEGG))#以数据框的形式展示KEGG的结果
browseKEGG(KEGG,"hsa04060")#打开KEGG官网查看通路的详细信息



KEGG_results <- as.data.frame(KEGG)######利用这段代码就可以获得我们想要的P<0.05的代码啦

####################单独想提取出我想要的那几行
KEGG_pathways<- KEGG_results[c(4,11,54,2,3,5,6,7,8,9),]   # 筛选我想要的通路
KEGG_pathways=arrange(KEGG_pathways,GeneRatio)##arrange函数的使用
KEGG_pathways$logPvalue <- (-log(KEGG_pathways$pvalue))



# 美化版本 --------------------------------------------------------------------
colnames(KEGG_pathways)

ggplot(KEGG_pathways, aes(x = Description, y = GeneRatio,fill = pvalue)) + #指定数据为merge_data，以Description作为x轴，以logPvalue作为y轴，以ONTOLOGY作为颜色填充。
  geom_bar(stat = "identity", width = 0.8, color = "black") + #添加一个柱状图层，使用identity统计方法（不进行任何变换），每个柱子的宽度为0.8，颜色为黑色
  coord_flip() + #翻转坐标轴，交换x轴和y轴，使得柱状图变为横向。
  #scale_y_continuous(expand = c(0,0))+ #设置y轴的范围为连续型变量（continuous），设置expand参数为c(0,0)，意味着y轴不需要扩展（不需要增加额外的空白）。
  theme(panel.background = element_rect(fill = "white"))+ #设置背景为白色，使用element_rect()函数设置元素的矩形形状
  labs(x = "KEGG pathways", y = "-log P-value")



#########################################GO_KEGG富集圈
library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(GOplot)#可视化
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")

# 导入差异分析的结果 ---------------------------------------------------------------
DEG <- read.csv("1.2DEG_deseq2(mRNAlog2(1.5)).csv",row.names = 1)
DEG=DEG[rownames(DEG)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]###这里我用全部的32个基因来做【因为富集分析也可以用指定的基因集来做】

DEG_use <- DEG

colnames(DEG_use) # 查看列名，下方选择列的依据
colnames(DEG_use)[1] <- "SYMBOL"# 修改列名，为下一步合并做准备
DEG_use <- DEG_use[,c("SYMBOL","log2FoldChange")] #只要基因名称和差异倍数两列




# 读取数据：基因差异分析的结果 ----------------------------------------------------------
Genes <- bitr(DEG_use$SYMBOL,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

data_use <- merge(DEG_use,Genes,by = "SYMBOL") #按照SYMBOL列合并数据，删除没有匹配到ENTREZID的基因行

# GO富集分析与整理 -------------------------------------------------------------
GO <- enrichGO(gene = data_use$ENTREZID, # 输入基因的 ENTREZID
               OrgDb = org.Hs.eg.db, # 注释信息来自 org.Hs.eg.db 数据库
               keyType = "ENTREZID", # 指定使用的基因标识符类型
               ont = "all", # 本次富集分析所使用的本体类型: "BP" (生物学过程)、"CC" (细胞组分)、"MF" (分子功能) 或 "all"（所有三个）
               pAdjustMethod = "BH", # 用于多重检验校正的 p 值校正方法 (Benjamini-Hochberg 方法)
               pvalueCutoff = 1, # p 值的阈值 (p 值低于此阈值的基因将被视为显著)
               qvalueCutoff = 1, # q 值的阈值 (q 值低于此阈值的基因将被视为显著)
               minGSSize = 5, # 基因集的最小大小 (基因的数量)
               maxGSSize = 5000, # 基因集的最大大小
               readable = TRUE) # 是否将 ENTREZID 转换为基因符号以提高可读性


# 将 GO 结果转换为数据框格式
GO_result <- as.data.frame(GO) 
# 选择显著富集功能通路，根据给定的阈值（Threshold）
GO_result <- GO_result[(GO_result$pvalue < 0.05 & GO_result$p.adjust < 0.05),] 
# 创建一个数据框，包含富集分析结果的GO类型
go_result <- data.frame(Category = GO_result$ONTOLOGY, 
                        ID = GO_result$ID, # 包含富集分析结果的通路 ID
                        Term = GO_result$Description, # 包含富集分析结果的通路描述
                        Genes = gsub("/", ", ", GO_result$geneID), # 包含富集分析结果中的基因 ID，多个 ID 用逗号分隔
                        adj_pval = GO_result$p.adjust) # 包含富集分析结果中校正后的 p 值



#########################################################画一个KEGG图
# 导入差异分析的结果 ---------------------------------------------------------------
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
DEG <- read.csv("1.2DEG_deseq2(mRNAlog2(1.5)).csv",row.names = 1)
DEG=DEG[rownames(DEG)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]###这里我用全部的32个基因来做【因为富集分析也可以用指定的基因集来做】

DEG_use <- DEG

colnames(DEG_use) # 查看列名，下方选择列的依据
colnames(DEG_use)[1] <- "SYMBOL"# 修改列名，为下一步合并做准备
DEG_use <- DEG_use[,c("SYMBOL","log2FoldChange")] #只要基因名称和差异倍数两列


# 读取数据：基因差异分析的结果 ----------------------------------------------------------
Genes <- bitr(DEG_use$SYMBOL,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

data_use <- merge(DEG_use,Genes,by = "SYMBOL") #按照SYMBOL列合并数据，删除没有匹配到ENTREZID的基因行

# KEGG富集分析与整理 -------------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# 将 KEGG 结果转换为数据框格式
KEGG_results <- DOSE::setReadable(KEGG,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENTREZID") #将ENTREZID转换为genesymbol
KEGG_result=as.data.frame(KEGG_results)
# 选择显著富集功能通路，根据给定的阈值（Threshold）
KEGG_result <- KEGG_result[(KEGG_result$pvalue < 0.05 & KEGG_result$p.adjust < 0.05),] 
# 创建一个数据框，包含富集分析结果的GO类型
kegg_result <- data.frame(ID = KEGG_result$ID, # 包含富集分析结果的通路 ID
                          Term = KEGG_result$Description, # 包含富集分析结果的通路描述
                          Genes = gsub("/", ", ", KEGG_result$geneID), # 包含富集分析结果中的基因 ID，多个 ID 用逗号分隔
                          adj_pval = KEGG_result$p.adjust) # 包含富集分析结果中校正后的 p 值




setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
R1=read.table("50.2hsa-miR-339-3p.txt")
R1=R1[-1,]
R2=read.table("50.3hsa-miR-17-5p.txt")
R2=R2[-1,]
R3=read.table("50.4hsa-miR-20b-5p.txt")
R3=R3[-1,]
R4=read.table("50.5 hsa-miR-210-3p.txt")
R4=R4[-1,]
R5=read.table("50.6hsa-miR-3187-3p.txt")
R5=R5[-1,]
R6=read.table("50.7hsa-miR-129-2-3p.txt")
R6=R6[-1,]
R7=read.table("50.8hsa-miR-423-5p.txt")
R7=R7[-1,]
R8=read.table("50.9hsa-miR-93-5p.txt")
R8=R8[-1,]
R9=read.table("50.10hsa-miR-9-3p.txt")
R9=R9[-1,]
R10=read.table("50.11hsa-miR-877-5p.txt")
R10=R10[-1,]
R11=read.table("50.12hsa-miR-887-3p.txt")
R11=R11[-1,]
R12=read.table("50.13hsa-miR-941.txt")
R12=R12[-1,]
R13=read.table("50.2hsa-miR-339-3p.txt")
R13=R13[-1,]

##################################合并
R=rbind(R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13)
R=R[!duplicated(R$V4),]
R1=R[R$V24!=0,]
write.csv(R1,"50.14去重复后预测的mRNA.csv")


####################################提取出mRNA
setwd("D://R/TCGA下载+整理数据/COAD数据/")
exp=read.csv("output_mRNA_lncRNA_expr/TCGA-COAD_mrna_expr_fpkm.csv",row.names = 1)
#### 过滤标准1：仅去除在所有样本里表达量都为零的基因（这是最低的过率标准）
exp = exp[rowSums(exp)>0,]##都为0
nrow(exp)
############################fpkm表达矩阵的处理
table(substr(colnames(exp),14,16))###查看肿瘤和正常样本的数目
Tumor <- grep('01A',colnames(exp))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_exp=exp[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
Normal <- grep('11A',colnames(exp))
Normal #正常样本所处位置
Normal_exp=exp[,Normal]
fpkm_exp=cbind(Normal_exp,Tumor_exp)
######提取出fpkm矩阵中的相关lncRNA
fpkm_exp1=fpkm_exp
############不要忘了log2的放缩
range(fpkm_exp1)
fpkm_exp1=log2(fpkm_exp1+1)###这里已经放缩了

###################临床信息与表达矩阵的结合
#表达矩阵
table(substr(colnames(fpkm_exp1),14,16))#查看肿瘤样本和正常样本数量
##只保留肿瘤样本，把正常的样本剔除掉
Tumor <- grep('01A',colnames(fpkm_exp1))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_lnc=fpkm_exp1[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
library(stringr)
expr_use = Tumor_lnc[,sort(colnames(Tumor_lnc))]
k = !duplicated(str_sub(colnames(expr_use),1,12))
table(k)
expr_use = expr_use[,k] 
colnames(expr_use) <- substr(colnames(expr_use),1,12)
##临床信息
setwd("D://R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
survival_data=read.csv("1.3结肠癌临床数据.csv",row.names = 1)
survival_data=survival_data[c(1,5,2,3,4,6,7,8,9)]
table(substr(rownames(survival_data),14,16))#查看肿瘤样本和正常样本数量
######只要肿瘤数据
Tumor <- grep('01A',rownames(survival_data))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
survival_data=survival_data[Tumor,]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
survival_data = survival_data[sort(rownames(survival_data)),]
library(stringr)
k = !duplicated(str_sub(rownames(survival_data),1,12))
table(k)
survival_data = survival_data[k,] 
rownames(survival_data) <- substr(rownames(survival_data),1,12)
library(stringr)
rownames(survival_data)=str_replace_all(rownames(survival_data),"-",".")
#将表达矩阵和生存信息合并
datExpr <- t(expr_use)
datExpr=as.data.frame(datExpr)
dim(datExpr)
dim(survival_data)
head(rownames(datExpr))
head(rownames(survival_data))
s = intersect(rownames(datExpr),rownames(survival_data));length(s)
####通过取交集，进行匹配【有的数据是只有表达数据，却没有临床信息，这样的样本是没有意义的】
datExpr1 = datExpr[s,]
survival_data2 = survival_data[s,]
dim(datExpr1)
dim(survival_data2)
identical(rownames(survival_data2),rownames(datExpr1))
#生存信息与表达矩阵合并
datExpr <- cbind(survival_data2,datExpr1)####只有这样才能够保持time和event在最前面
#保存数据
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
write.csv(datExpr,"50.15 临床数据与全部的mRNAs表达数据结合.csv")

######################################确保lncRNA与mRNA成正比
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
exp1=read.csv("50.15 临床数据与全部的mRNAs表达数据结合.csv",row.names = 1)
exp1=exp1[,-c(1:9)]
exp2=read.csv("44.1 临床数据与4个预后基因lncRNAs表达数据结合.csv",row.names = 1)
exp2=exp2[,-c(1:9)]

########################################矩阵的合并
exp5=cbind(exp2,exp1)
range(exp5)####计算相关性可以暂时不用log2缩小差距

#################################批量计算基因间的相关性
# 计算相关性R和显著性p值 ------------------------------------------------------------
data_use=exp5
#"LINC01133","AP003555.1"
target_gene <- c("LINC01133")###只有这两个预测到了mRNA

#################################批量计算基因间的相关性
#生成一个空白的数据框，存放结果
result <- c()
#创建要进行相关性分析的基因列表
gene_list <- colnames(data_use)[1:19490]#进行相关性分析的范围【我们是选择全部的基因】
#循环读取
for (gene in gene_list) {
  for(target in target_gene){
    message(gene,target)#打印【报出信息】
    cor_R <- cor.test(x = data_use[,target],y = data_use[,gene],method = 'pearson')###基因近似正态分布
    P.value<-cor_R[["p.value"]]
    cor=cor_R[["estimate"]][["cor"]]
    temp_result<-c(gene,target,cor,P.value)#让这个变量一直在变
    result <- rbind(result, temp_result)
  }
  
}
colnames(result)=c("mRNA","LINC01133","cor_R","cor_P")
result=as.data.frame(result)
result=na.omit(result)###会出现缺失值
result$cor_R=as.numeric(result$cor_R)
result$cor_P=as.numeric(result$cor_P)
#自行设置阈值（二选一），筛选要的基因对
table(result$cor_R >0.3)#查看相关性R大于0.3的基因个数
table(result$cor_P < 0.05)#查看p值小于0.05（阈值自选）的基因个数
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
write.csv(result,"50.16LINC01133_mRNA的分析.csv")

index <- result$cor_P < 0.05
index1<-result$cor_R >0.3#注意：正常情况——相关性系数在【-1，1】之间都是相关性比较高的，但是我们这里选择的是与基因正相关的
result_filt <- result[index&index1,]
########去重复
library(dplyr)
result_filt1 <- result_filt[!duplicated(result_filt$mRNA),]
rownames(result_filt1)=result_filt1$mRNA
write.csv(result_filt1,"50.16LINC01133_mRNA的分析(0.3).csv")




###########################################################microRNA
setwd("D://R/TCGA下载+整理数据/COAD数据/output_miRNA_expr/")
load("TCGA-COAD_mirna_expr_rpm.rdata")
A=as.data.frame(rownames(mirna_expr_rpm))
###################
exp=mirna_expr_rpm[rownames(mirna_expr_rpm)%in%c("hsa-mir-210","hsa-mir-3187","hsa-mir-129-2",
                                                 "hsa-mir-941-1","hsa-mir-941-2","hsa-mir-941-3",
                                                 "hsa-mir-941-4","hsa-mir-941-5",
                                                 "hsa-mir-4435-1","hsa-mir-4435-2","hsa-mir-4745",
                                                 "hsa-mir-887","hsa-mir-128-1","hsa-mir-15a",
                                                 "hsa-mir-17","hsa-mir-186","hsa-mir-20b",
                                                 "hsa-mir-339","hsa-mir-423","hsa-mir-877",
                                                 "hsa-mir-9-1","hsa-mir-9-2","hsa-mir-9-3","hsa-mir-93"),]

############################fpkm表达矩阵的处理
table(substr(colnames(exp),14,16))###查看肿瘤和正常样本的数目
Tumor <- grep('01A',colnames(exp))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_exp=exp[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
Normal <- grep('11A',colnames(exp))
Normal #正常样本所处位置
Normal_exp=exp[,Normal]
fpkm_exp=cbind(Normal_exp,Tumor_exp)
######提取出fpkm矩阵中的相关lncRNA
fpkm_exp1=fpkm_exp
############不要忘了log2的放缩
range(fpkm_exp1)
fpkm_exp1=log2(fpkm_exp1+1)###这里已经放缩了

###################临床信息与表达矩阵的结合
#表达矩阵
table(substr(colnames(fpkm_exp1),14,16))#查看肿瘤样本和正常样本数量
##只保留肿瘤样本，把正常的样本剔除掉
Tumor <- grep('01A',colnames(fpkm_exp1))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_lnc=fpkm_exp1[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
library(stringr)
expr_use = Tumor_lnc[,sort(colnames(Tumor_lnc))]
k = !duplicated(str_sub(colnames(expr_use),1,12))
table(k)
expr_use = expr_use[,k] 
colnames(expr_use) <- substr(colnames(expr_use),1,12)
##临床信息
setwd("D://R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
survival_data=read.csv("1.3结肠癌临床数据.csv",row.names = 1)
survival_data=survival_data[c(1,5,2,3,4,6,7,8,9)]
table(substr(rownames(survival_data),14,16))#查看肿瘤样本和正常样本数量
######只要肿瘤数据
Tumor <- grep('01A',rownames(survival_data))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
survival_data=survival_data[Tumor,]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
survival_data = survival_data[sort(rownames(survival_data)),]
library(stringr)
k = !duplicated(str_sub(rownames(survival_data),1,12))
table(k)
survival_data = survival_data[k,] 
rownames(survival_data) <- substr(rownames(survival_data),1,12)
library(stringr)
rownames(survival_data)=str_replace_all(rownames(survival_data),"-",".")
#将表达矩阵和生存信息合并
datExpr <- t(expr_use)
datExpr=as.data.frame(datExpr)
dim(datExpr)
dim(survival_data)
head(rownames(datExpr))
head(rownames(survival_data))
rownames(datExpr)=str_replace_all(rownames(datExpr),"-",".")
s = intersect(rownames(datExpr),rownames(survival_data));length(s)
####通过取交集，进行匹配【有的数据是只有表达数据，却没有临床信息，这样的样本是没有意义的】
datExpr1 = datExpr[s,]
survival_data2 = survival_data[s,]
dim(datExpr1)
dim(survival_data2)
identical(rownames(survival_data2),rownames(datExpr1))
#生存信息与表达矩阵合并
datExpr <- cbind(survival_data2,datExpr1)####只有这样才能够保持time和event在最前面
#保存数据
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
write.csv(datExpr,"51.1 临床数据与预测的microRNAs表达数据结合.csv")

######################################确保lncRNA与mRNA成正比
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
exp1=read.csv("51.11 临床数据与预测的microRNAs表达数据结合.csv",row.names = 1)
exp1=exp1[,-c(1:9)]
exp2=read.csv("50.15 临床数据与预测的mRNAs表达数据结合.csv",row.names = 1)
exp2=exp2[,-c(1:9)]
exp3=read.csv("50.18最终预测的mRNA.csv",row.names = 1)
exp2=exp2[,exp3$mRNA]
exp4=read.csv("50.14去重复后预测的mRNA.csv",row.names = 1)
########################################矩阵的合并
exp2=exp2[rownames(exp1),]
exp5=cbind(exp2,exp1)
range(exp5)####计算相关性可以暂时不用log2缩小差距

#################################批量计算基因间的相关性
# 计算相关性R和显著性p值 ------------------------------------------------------------
data_use=exp5
colnames(exp1)
#"hsa.mir.129.2" "hsa.mir.17"    "hsa.mir.20b"   "hsa.mir.210"   "hsa.mir.3187" 
#"hsa.mir.339"   "hsa.mir.423"   "hsa.mir.877"   "hsa.mir.887"   "hsa.mir.9.3"  
#"hsa.mir.93"   

target_gene <- c("hsa.mir.93")

#################################批量计算基因间的相关性
#生成一个空白的数据框，存放结果
result <- c()
#创建要进行相关性分析的基因列表
gene_list <- colnames(data_use)[1:84]#进行相关性分析的范围【我们是选择全部的基因】
#循环读取
for (gene in gene_list) {
  for(target in target_gene){
    message(gene,target)#打印【报出信息】
    cor_R <- cor.test(x = data_use[,target],y = data_use[,gene],method = 'pearson')###基因近似正态分布
    P.value<-cor_R[["p.value"]]
    cor=cor_R[["estimate"]][["cor"]]
    temp_result<-c(gene,target,cor,P.value)#让这个变量一直在变
    result <- rbind(result, temp_result)
  }
  
}
colnames(result)=c("mRNA","microRNA","cor_R","cor_P")
result=as.data.frame(result)
result=na.omit(result)###会出现缺失值
result$cor_R=as.numeric(result$cor_R)
result$cor_P=as.numeric(result$cor_P)
rownames(result)=result$mRNA


exp5=exp4[exp4$V2=="hsa-miR-93-5p",]
K=intersect(exp5$V4,result$mRNA);K
result1=result[K,]
table(result1$cor_R <0)
table(result1$cor_P < 0.05)
index <- result1$cor_P < 0.05
index1<-result1$cor_R <0
result1 <- result1[index&index1,]
write.csv(result1,"51.7 hsa-miR-423-5p_mRNA_microRNA.csv")

###############################################提取出最终的microRNA
M1=read.csv("51.3 hsa-miR-129-2-3p_mRNA_microRNA.csv",row.names = 1)
M2=read.csv("51.4 hsa-miR-17-5p_mRNA_microRNA.csv",row.names = 1)
M3=read.csv("51.5 hsa-miR-3187-3p_mRNA_microRNA.csv",row.names = 1)
M4=read.csv("51.6 hsa-miR-339-3p_mRNA_microRNA.csv",row.names = 1)
M5=read.csv("51.7 hsa-miR-423-5p_mRNA_microRNA.csv",row.names = 1)
M=rbind(M1,M2,M3,M4,M5)
M1=M[M$cor_R< -0.3,]
write.csv(M,"51.8 全部的microRNA.csv")
write.csv(M1,"51.9 小于0.3的microRNA.csv")


#############################################lncRNA_microRNA_mRNA结合
L=read.csv("50.18最终预测的mRNA.csv",row.names = 1)
M=read.csv("51.9 小于0.3的microRNA.csv",row.names = 1)
rownames(L)=L$mRNA
rownames(M)=M$mRNA
K=intersect(L$mRNA,M$mRNA)
L1=L[K,]
M1=M[K,]
L1$microRNA=M1$microRNA
write.csv(L1,"51.10 最终的ceRNA.csv")

###########################################再做一个富集分析
library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(ggplot2)#分面绘图所需
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
# 读取数据：基因差异分析的结果 ----------------------------------------------------------
upGene <- read.csv("1.2DEG_deseq2(mRNAlog2(1.5)).csv",row.names = 1)
table(upGene$regulate)
L1=read.csv("51.10 最终的ceRNA.csv",row.names = 1)
#upGene=upGene[upGene$regulate!="unchanged",]#GSEA富集分析不需要寻找差异基因
upGene2=upGene[rownames(L1),]
#upGene3=upGene2[upGene2$regulate!="unchanged",]
upGene=upGene2
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

# GO富集分析与结果存储【GO没有结果】 -------------------------------------------------------------
GO <- enrichGO(gene = Genes$ENTREZID, #输入基因的"ENTREZID"
               OrgDb = org.Hs.eg.db,#注释信息
               keyType = "ENTREZID",
               ont = "all",     #可选条目BP/CC/MF
               pAdjustMethod = "BH", #p值的校正方式
               pvalueCutoff = 0.05,   #pvalue的阈值
               qvalueCutoff = 0.05, #qvalue的阈值
               minGSSize = 5,
               maxGSSize = 5000,
               readable = TRUE)   #是否将entrez id转换为symbol
AA=GO@result
AA=AA[AA$pvalue<0.05,]
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
write.csv(AA,"52.1 预测的mRNAGO富集分析.csv")


# KEGG富集分析与结果存储【没有结果】 -----------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)

####结果存储
KEGG_results <- DOSE::setReadable(KEGG,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENTREZID") #将ENTREZID转换为genesymbol
AA=KEGG_results@result
AA=AA[AA$pvalue<0.05,]
write.csv(AA,"52.2 预测的mRNA的KEGG富集分析.csv")
miRNA1=read.table("50.2 hsa-miR-17-5p.txt",header = T,skip = 4)
miRNA2=read.table("50.3 hsa-miR-339-3p.txt",header = T,skip = 4)
miRNA3=read.table("50.4 hsa-miR-93-5p.txt",header = T,skip = 4)
miRNA4=read.table("50.5 hsa-miR-20b-5p.txt",header = T,skip = 4)
miRNA5=read.table("50.6 hsa-miR-423-5p.txt",header = T,skip = 4)
miRNA6=read.table("50.7 hsa-miR-877-5p.txt",header = T,skip = 4)
miRNA7=read.table("50.8 hsa-miR-9-3p.txt",header = T,skip = 4)
miRNA8=read.table("50.9 hsa-miR-129-2-3p.txt",header = T,skip = 4)
miRNA9=read.table("50.10 hsa-miR-210-3p.txt",header = T,skip = 4)
miRNA10=read.table("50.11 hsa-miR-3187-3p.txt",header = T,skip = 4)
miRNA11=read.table("50.12 hsa-miR-887-3p.txt",header = T,skip = 4)
miRNA12=read.table("50.13 hsa-miR-941.txt",header = T,skip = 4)
miRNA=rbind(miRNA1,miRNA2,miRNA3,miRNA4,miRNA5,miRNA6,miRNA7,
            miRNA8,miRNA9,miRNA10,miRNA11,miRNA12)
miRNA=miRNA[!duplicated(miRNA$geneName),]
write.csv(miRNA,"54.5 所有miRNA整合后去重复.csv")
miRNA1=read.table("50.13 hsa-miR-941.txt",header = T,skip = 4)
miRNA1=miRNA1[!duplicated(miRNA1$geneName),]
write.table(miRNA,"50.13 hsa-miR-941.txt",sep = "/t",)
?write.table
library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(ggplot2)#分面绘图所需
library(ggsci)
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
# 读取数据：基因差异分析的结果 ----------------------------------------------------------
upGene <- read.csv("1.2DEG_deseq2(mRNAlog2(1.5)).csv",row.names = 1)
upGene$regulate <- ifelse(upGene$padj > 0.05, "unchanged",
                          ifelse(upGene$log2FoldChange > 1, "up-regulated",
                                 ifelse(upGene$log2FoldChange < -1, "down-regulated", "unchanged")))
table(upGene$regulate)
upGene2=upGene[rownames(upGene)%in%c(
  "PRRX1","ITPRIP","TRIM8","AHNAK","ETS1","CHST11","FICD",
  "OAS2","CX3CL1","CCDC88A","CYTIP","EPHA4","NXPE3","CD109",
  "ESR1","SAMD9L","CAV2"),]###用我们预测的17个基因来做
upGene=upGene2[upGene2$regulate!="unchanged",]
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

# GO富集分析与结果存储【GO没有结果】 -------------------------------------------------------------
GO <- enrichGO(gene = Genes$ENTREZID, #输入基因的"ENTREZID"
               OrgDb = org.Hs.eg.db,#注释信息
               keyType = "ENTREZID",
               ont = "all",     #可选条目BP/CC/MF
               pAdjustMethod = "BH", #p值的校正方式
               pvalueCutoff = 0.05,   #pvalue的阈值
               qvalueCutoff = 0.05, #qvalue的阈值
               minGSSize = 5,
               maxGSSize = 5000,
               readable = TRUE)   #是否将entrez id转换为symbol
AA=GO@result
AA=AA[AA$pvalue<0.05,]
write.csv(AA,"55.5 lncRNA_microRNA_mRNA的GO富集分析.csv")


# 或许你想展示某些特定条目 ------------------------------------------------------------
GO_result <- as.data.frame(GO)

####################单独想提取出我想要的那几行
BP_top10<- GO_result[c(1,4,70,58,62),]   # 筛选 ONTOLOGY 为 BP 的子集
BP_top10%>%arrange(BP_top10$pvalue) 
BP_top10$logPvalue <- (-log(BP_top10$pvalue))####一定要转化为这个样子，这样才可以画图更好看
# 美化版本 --------------------------------------------------------------------
colnames(BP_top10)
BP_top10$Description <- factor(BP_top10$Description)
pdf("GO.pdf",width = 18,height = 16)
ggplot(BP_top10, aes(x = Description, y = logPvalue,fill=Description)) + #指定数据为merge_data，以Description作为x轴，以logPvalue作为y轴，以ONTOLOGY作为颜色填充。
  geom_bar(stat = "identity", width = 0.8, color = "black") + #添加一个柱状图层，使用identity统计方法（不进行任何变换），每个柱子的宽度为0.8，颜色为黑色
  scale_x_discrete(limits = unique(BP_top10$Description)) + #设置x轴的标签顺序，根据ONTOLOGY对Description排序，以便按照ONTOLOGY的顺序显示。unique()函数保证每个标签只出现一次。
  coord_flip() + #翻转坐标轴，交换x轴和y轴，使得柱状图变为横向。
  scale_y_continuous(expand = c(0,0))+ #设置y轴的范围为连续型变量（continuous），设置expand参数为c(0,0)，意味着y轴不需要扩展（不需要增加额外的空白）。
  theme(panel.background = element_rect(fill = "white"))+ #设置背景为白色，使用element_rect()函数设置元素的矩形形状
  theme(axis.text = element_text(size = 16, family = "Arial"),
        axis.title = element_text(size=18,family="Arial",colour = "black"),
        legend.text = element_text(size = 16, family = "Arial"),
        legend.title = element_text(size = 18, family = "Arial"),
        legend.position = "top")+ #设置各种文本元素的字体、大小和颜色，包括坐标轴文本、坐标轴标题、图例文本和图例标题。将图例放置在顶部。
  labs(x = "GO term", y = "-log P-value")+
  scale_fill_nejm(alpha = 0.8) #使用NEJM（New England Journal of Medicine）颜色调色板，设置颜色透明度为0.8，对填充颜色进行比例缩放。


# KEGG富集分析与结果存储【没有结果】 -----------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

####结果存储
KEGG_results <- DOSE::setReadable(KEGG,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENTREZID") #将ENTREZID转换为genesymbol
KEGG_results <- as.data.frame(KEGG)
write.csv(KEGG_results,"43.2 凋亡基因的KEGG富集分析.csv")

#################################################基础绘图
barplot(KEGG,
        x = "GeneRatio",
        color = "pvalue",
        showCategory = 15,
        title = "KEGG_enrichment")#标题
dotplot(KEGG)
View(as.data.frame(KEGG))#以数据框的形式展示KEGG的结果
browseKEGG(KEGG,"hsa04060")#打开KEGG官网查看通路的详细信息



KEGG_results <- as.data.frame(KEGG)######利用这段代码就可以获得我们想要的P<0.05的代码啦

####################单独想提取出我想要的那几行
KEGG_pathways<- KEGG_results[c(4,11,54,2,3,5,6,7,8,9),]   # 筛选我想要的通路
KEGG_pathways=arrange(KEGG_pathways,GeneRatio)##arrange函数的使用
KEGG_pathways$logPvalue <- (-log(KEGG_pathways$pvalue))



# 美化版本 --------------------------------------------------------------------
colnames(KEGG_pathways)

ggplot(KEGG_pathways, aes(x = Description, y = GeneRatio,fill = pvalue)) + #指定数据为merge_data，以Description作为x轴，以logPvalue作为y轴，以ONTOLOGY作为颜色填充。
  geom_bar(stat = "identity", width = 0.8, color = "black") + #添加一个柱状图层，使用identity统计方法（不进行任何变换），每个柱子的宽度为0.8，颜色为黑色
  coord_flip() + #翻转坐标轴，交换x轴和y轴，使得柱状图变为横向。
  #scale_y_continuous(expand = c(0,0))+ #设置y轴的范围为连续型变量（continuous），设置expand参数为c(0,0)，意味着y轴不需要扩展（不需要增加额外的空白）。
  theme(panel.background = element_rect(fill = "white"))+ #设置背景为白色，使用element_rect()函数设置元素的矩形形状
  labs(x = "KEGG pathways", y = "-log P-value")



#########################################GO_KEGG富集圈
library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(GOplot)#可视化
setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")

# 导入差异分析的结果 ---------------------------------------------------------------
DEG <- read.csv("1.2DEG_deseq2(mRNAlog2(1.5)).csv",row.names = 1)
DEG=DEG[rownames(DEG)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]###这里我用全部的32个基因来做【因为富集分析也可以用指定的基因集来做】

DEG_use <- DEG

colnames(DEG_use) # 查看列名，下方选择列的依据
colnames(DEG_use)[1] <- "SYMBOL"# 修改列名，为下一步合并做准备
DEG_use <- DEG_use[,c("SYMBOL","log2FoldChange")] #只要基因名称和差异倍数两列




# 读取数据：基因差异分析的结果 ----------------------------------------------------------
Genes <- bitr(DEG_use$SYMBOL,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

data_use <- merge(DEG_use,Genes,by = "SYMBOL") #按照SYMBOL列合并数据，删除没有匹配到ENTREZID的基因行

# GO富集分析与整理 -------------------------------------------------------------
GO <- enrichGO(gene = data_use$ENTREZID, # 输入基因的 ENTREZID
               OrgDb = org.Hs.eg.db, # 注释信息来自 org.Hs.eg.db 数据库
               keyType = "ENTREZID", # 指定使用的基因标识符类型
               ont = "all", # 本次富集分析所使用的本体类型: "BP" (生物学过程)、"CC" (细胞组分)、"MF" (分子功能) 或 "all"（所有三个）
               pAdjustMethod = "BH", # 用于多重检验校正的 p 值校正方法 (Benjamini-Hochberg 方法)
               pvalueCutoff = 1, # p 值的阈值 (p 值低于此阈值的基因将被视为显著)
               qvalueCutoff = 1, # q 值的阈值 (q 值低于此阈值的基因将被视为显著)
               minGSSize = 5, # 基因集的最小大小 (基因的数量)
               maxGSSize = 5000, # 基因集的最大大小
               readable = TRUE) # 是否将 ENTREZID 转换为基因符号以提高可读

setwd("D:/R/SLC12A9——结肠癌/泛凋亡弄出来了！！/")
DEG <- read.csv("1.2DEG_deseq2(mRNAlog2(1.5)).csv",row.names = 1)
DEG=DEG[rownames(DEG)%in%c(
  "PYCARD","CASP3","CASP8","RIPK1","FADD","GSDME","CASP7",
  "CASP1","CASP11","TAK1","MLKL","GSDMD","NLRP3","ZBP1",
  "TNF","IFNG","ADAR","IL1B","NLRC4","IRF8","CASP4","CASP5",
  "IL18","NLRP12","DDX3X","CASP6","PARP1","NAIP","AIM2","IRF1",
  "TLR9","RIPK3","TRADD","NFS1"),]###这里我用全部的32个基因来做【因为富集分析也可以用指定的基因集来做】

DEG_use <- DEG

colnames(DEG_use) # 查看列名，下方选择列的依据
colnames(DEG_use)[1] <- "SYMBOL"# 修改列名，为下一步合并做准备
DEG_use <- DEG_use[,c("SYMBOL","log2FoldChange")] #只要基因名称和差异倍数两列


# 读取数据：基因差异分析的结果 ----------------------------------------------------------
Genes <- bitr(DEG_use$SYMBOL,  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

data_use <- merge(DEG_use,Genes,by = "SYMBOL") #按照SYMBOL列合并数据，删除没有匹配到ENTREZID的基因行

# KEGG富集分析与整理 -------------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)





exp1=read.csv("TCGA-COAD_mrna_expr_fpkm.csv",row.names = 1)
exp2=read.csv("TCGA-COAD_lncrna_expr_fpkm.csv",row.names = 1)
exp=rbind(exp1,exp2)
table(substr(colnames(exp),14,16))
Tumor <- grep('01A',colnames(exp))
Tumor 
Tumor_exp=exp[,Tumor]
Normal <- grep('11A',colnames(exp))
Normal #正常样本所处位置
Normal_exp=exp[,Normal]
fpkm_exp=cbind(Normal_exp,Tumor_exp)#####重点：这里要求肿瘤样本和正常组的样本都要有的
Expr=t(apply(fpkm_exp,1,function(x){x-(mean(x))}))###这样标准化后的值有正有负

exp=read.csv("56.3 TIDE结果.csv",row.names = 1)
exp=t(exp)
exp=as.data.frame(exp)
table(substr(colnames(exp),14,16))###查看肿瘤和正常样本的数目
Tumor <- grep('01A',colnames(exp))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_exp=exp[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
Normal <- grep('11A',colnames(exp))
Normal #正常样本所处位置
Normal_exp=exp[,Normal]
fpkm_exp=cbind(Normal_exp,Tumor_exp)
######提取出fpkm矩阵中的相关lncRNA
fpkm_exp1=fpkm_exp
###################临床信息与表达矩阵的结合
#表达矩阵
table(substr(colnames(fpkm_exp1),14,16))#查看肿瘤样本和正常样本数量
##只保留肿瘤样本，把正常的样本剔除掉
Tumor <- grep('01A',colnames(fpkm_exp1))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
Tumor_lnc=fpkm_exp1[,Tumor]###为了数据的严谨性，我们只需要01A的癌症组织
#表达矩阵以样本为核心去重
library(stringr)
expr_use = Tumor_lnc[,sort(colnames(Tumor_lnc))]
k = !duplicated(str_sub(colnames(expr_use),1,12))
table(k)
expr_use = expr_use[,k] 
colnames(expr_use) <- substr(colnames(expr_use),1,12)
survival_data=read.csv("1.3结肠癌临床数据.csv",row.names = 1)
survival_data=survival_data[c(1,5,2,3,4,6,7,8,9)]
table(substr(rownames(survival_data),14,16))#查看肿瘤样本和正常样本数量
Tumor <- grep('01A',rownames(survival_data))
Tumor  #肿瘤样本所处位置
#01A：癌症组织   01B：福尔马林浸泡的组织 01C：其他 02A：复发组织 06A：转移组织
survival_data=survival_data[Tumor,]###为了数据的严谨性，我们只需要01A的癌症组织
survival_data = survival_data[sort(rownames(survival_data)),]
library(stringr)
k = !duplicated(str_sub(rownames(survival_data),1,12))
table(k)
survival_data = survival_data[k,] 
rownames(survival_data) <- substr(rownames(survival_data),1,12)
library(stringr)
rownames(survival_data)=str_replace_all(rownames(survival_data),"-",".")
datExpr <- t(expr_use)
datExpr=as.data.frame(datExpr)
dim(datExpr)
dim(survival_data)
head(rownames(datExpr))
head(rownames(survival_data))
s = intersect(rownames(datExpr),rownames(survival_data));length(s)
datExpr1 = datExpr[s,]
survival_data2 = survival_data[s,]
dim(datExpr1)
dim(survival_data2)
identical(rownames(survival_data2),rownames(datExpr1))
datExpr <- cbind(survival_data2,datExpr1)####只有这样才能够保持time和event在最前面




group=read.csv("29.2免疫的分组.csv",row.names = 1)
table(group$group)
datExpr=datExpr[rownames(group),]
identical(rownames(group),rownames(datExpr))
group$group=as.factor(group$group)
group$group=relevel(group$group,ref = "high")
group$TIDE=datExpr$TIDE
group$Dysfunction=datExpr$Dysfunction
group$Exclusion=datExpr$Exclusion

group1=group
group1$TIDE=as.numeric(group1$TIDE)
group1$Dysfunction=as.numeric(group1$Dysfunction)
group1$Exclusion=as.numeric(group1$Exclusion)

library(ggplot2)
library(ggpubr)

ggplot(data=group1, aes(x=group, y=Dysfunction,fill=group)) +
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  theme_bw()

