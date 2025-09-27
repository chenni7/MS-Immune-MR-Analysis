library(TwoSampleMR)
library(ggplot2)
library(foreach)


setwd("D://91血细胞//131_91种血细胞孟德尔随机化//bc91//1_step1")

iddf=read.table("bc91id.txt",header =T,sep = "\t")


bioid=as.vector(iddf$id)

result=data.frame()

foreach(i=bioid, .errorhandling = "pass") %do%{
 
  expo_rt<-read_exposure_data(
    filename = paste0("clump/",i,".txt"),
    sep = "\t",
    snp_col = "rsids2",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
    samplesize_col = "N")
  
  
  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "myoutcome.gz",
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval")
  
  harm_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=1)
  
  
  harm_rt<-subset(harm_rt,mr_keep==TRUE)
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2)) /
    (2 * (harm_rt$beta.exposure^2) +
       2 * harm_rt$samplesize.exposure * harm_rt$se.exposure^2)
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf<- mean( harm_rt$f)
  harm_rt<-harm_rt[harm_rt$f>10,]
  
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result) 
  if (mr_result$pval[3]<0.05){
    result=rbind(result,cbind(id=i,pvalue=result_or$pval[3])) 
    filename=paste0("result/",i)
    dir.create(filename) 
    write.table(harm_rt, file =paste0(filename,"/harmonise.txt"),row.names = F,sep = "\t",quote = F)
    write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
    pleiotropy=mr_pleiotropy_test(harm_rt)
    write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
    heterogeneity=mr_heterogeneity(harm_rt)
    write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
    
    
    p1 <- mr_scatter_plot(mr_result, harm_rt)
    ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
    
   
    
    singlesnp_res<- mr_singlesnp(harm_rt)
    singlesnpOR=generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
    p2 <- mr_forest_plot(singlesnp_res)
    ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
    sen_res<- mr_leaveoneout(harm_rt)
    p3 <- mr_leaveoneout_plot(sen_res)
    ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
    res_single <- mr_singlesnp(harm_rt)
    p4 <- mr_funnel_plot(singlesnp_res)
    ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
    presso=run_mr_presso(harm_rt,NbDistribution = 1000)
    capture.output(presso,file = paste0(filename,"/presso.txt"))
  }
}
write.table(result,"result.txt",sep = "\t",quote = F,row.names = F)







