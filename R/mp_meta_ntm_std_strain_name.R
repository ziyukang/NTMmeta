# MP meta analysis dataset preprocessing
# transform strain names to standard names

library(tidyverse)
library(readxl)
require(openxlsx)

ntm <- read_excel("MP meta分析数据收集表核查prev_eng.xlsx",sheet = "数据记录表", skip = 1)
ntm.name <- names(ntm)[23:119]

ntm.group <- as.data.frame(read_excel("MP meta分析数据收集表核查prev_eng.xlsx",sheet = "菌名校对和SGM_RGM分类"))

ntm.name.std <- ntm.group$菌株标准名称
names(ntm.name.std) <- ntm.group$菌名

ntm.std.data <- ntm[,1:22]
ntm.std.data$NumbervalidateNTM <- ntm$总计
for(name.std in ntm.name.std[!duplicated(ntm.name.std)]){
  name.data <- names(ntm.name.std[ntm.name.std == name.std])
  ntm.std.data[,name.std] <- rowSums(ntm[,name.data],na.rm=T)
}

write.xlsx(ntm.std.data, file = "MP meta分析数据收集表核查最终_标准菌株名.xlsx")

