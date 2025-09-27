library(TwoSampleMR)
library(ggplot2)
setwd("D://双细胞孟德尔随机化//133_双细胞孟德尔随机化//double_cell//4_median")

#######two step###########
#step 1 bc to outcome
harm_rt=read.table("harmonise1.txt",sep = "\t",header = T)

mr_result<- mr(harm_rt)
View(mr_result)
result_or=generate_odds_ratios(mr_result)
write.table(result_or[,5:ncol(result_or)],"OR1.txt",row.names = F,sep = "\t",quote = F)

#step 3 bc to imc
harm_rt3=read.table("harmonise2.txt",sep = "\t",header = T)
mr_result3<- mr(harm_rt3)
View(mr_result3)
result_or3=generate_odds_ratios(mr_result3)
write.table(result_or3[,5:ncol(result_or3)],"OR3.txt",row.names = F,sep = "\t",quote = F)

#step 4 imc to outcome
harm_rt4=read.table("harmonise3.txt",sep = "\t",header = T)
mr_result4<- mr(harm_rt4)
View(mr_result4)
result_or4=generate_odds_ratios(mr_result4)
write.table(result_or4[,5:ncol(result_or4)],"OR4.txt",row.names = F,sep = "\t",quote = F)

#####mediation effect
beta_all=mr_result[3,"b"]

beta1=mr_result3[3,"b"]

beta2=mr_result4[3,"b"]

#中介效应（间接效应）
beta12=beta1*beta2

#直接效应
beta_dir=beta_all-beta12

se=sqrt(mr_result3[3,"b"]^2*mr_result3[3,"se"]^2+mr_result4[3,"b"]^2*mr_result4[3,"se"]^2)

#Z
Z=beta12/se

#p
P=2*pnorm(q=abs(Z), lower.tail=FALSE)

#####95%可信区间
lci=beta12-1.96*se

uci=beta12+1.96*se

#中介效应所占的比例
beta12_p=beta12/beta_all
lci_p=lci/beta_all
uci_p=uci/beta_all

meresult=data.frame(beta_all,beta1,beta2,beta12,beta_dir,Z,P,lci,uci,beta12_p,lci_p,uci_p)
write.table(meresult,"meresult.txt",quote = F,sep = "\t",row.names = F)
