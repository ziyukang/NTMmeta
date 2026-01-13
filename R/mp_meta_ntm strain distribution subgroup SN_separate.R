# meta analysis of prevalence of NTM species 

library(tidyverse)
library(metafor) 
library(meta)
library(readxl)
library(r2rtf)
library(stringi) # 字符串处理
require(openxlsx)
require(DescTools)

set.seed(20251119)
format_pval <- function(pval){
  
  if(pval < 0.001) pval <- "< 0.001"
  else pval <- format(round(pval,3), nsmall = 3)
  return(pval)
}

ntm <- read_excel("MP meta分析数据收集表核查最终_标准菌株名.xlsx",skip = 0)
ntm.name <- names(ntm)[23:112]

#incorporate SGM and RGM information
ntm.group <- as.data.frame(read_excel("MP meta分析数据收集表核查prev_eng.xlsx",sheet = "菌名校对和SGM_RGM分类"))
rownames(ntm.group) <- ntm.group$菌名
ntm.group <- ntm.group[!duplicated(ntm.group$菌株标准名称),]
rownames(ntm.group) <- ntm.group$菌株标准名称

name.sgm <- ntm.group[ntm.group[,"RGM/SGM"] == "SGM","菌株标准名称"]
name.sgm <- name.sgm[!is.na(name.sgm)]
name.rgm <- ntm.group[ntm.group[,"RGM/SGM"] == "RGM","菌株标准名称"]
name.rgm <- name.rgm[!is.na(name.rgm)]

ntm$SGM <- rowSums(ntm[,name.sgm],na.rm=T)
ntm$RGM <- rowSums(ntm[,name.rgm],na.rm=T)

ntm <- ntm %>%
  select(FirstAuthor0, FirstAuthor, PublicationYear, StudyDesign, 
       StudyRegion, Regions, ScreeningPatientNumber, NumbervalidateNTM,
       prevalence, SGM, RGM, RegionSN) %>%
  mutate(Regions = factor(Regions, levels = c("东北","西北","中部","东南","西南"))) %>%
  filter(NumbervalidateNTM > 0) %>%
  filter(!is.na(NumbervalidateNTM)) %>%
  mutate(Study = paste(FirstAuthor, PublicationYear, sep = ",")) %>%
  mutate(ScreeningPatientNumber = as.numeric(ScreeningPatientNumber))

strain_res <- data.frame()
for(name in c("SGM","RGM")) {
  da <- ntm[, c("Study", "FirstAuthor", "PublicationYear", 
                "StudyRegion", "NumbervalidateNTM", name, "RegionSN")]
  if(sum(!da[,name] == 0) <= 0) next
  da <- da[!da[,name] == 0,]
  
  #Only choose one of the three transformation methods 
  ies.da <- escalc(xi = unlist(da[,name]), ni = da$NumbervalidateNTM, add = 0, measure = "PFT")
  ies.da$StudyRegion <- da$StudyRegion
  pes.da <- rma(yi, vi, data = ies.da, method = "REML")
  pes <- predict(pes.da, transf = transf.ipft.hm, targ = list(ni=1/(pes.da$se)^2))          
  
  num_study <- pes.da$k
  n <- sum(da[,name])
  N <- sum(da$NumbervalidateNTM)
  n_N <- paste0(n,"/", N)
  Q.pval <- format_pval(pes.da$QEp)
  I2 <- format(round(pes.da$I2,1), nsmall = 1)
  est <- format(round(100*c(pes[[1]], pes[[3]], pes[[4]]),1), nsmall = 1, trim = TRUE)
  strain.prev <- paste0(est[1], "(", est[2], "-", est[3], ")")

  if(nrow(da) == 1){
    Q.pval <- "-"
    I2 <- "-"
    est <- format(round(100*BinomCI(n,N,method="clopper-pearson"),1), nsmall = 1)
    strain.prev <- paste0(est[1], "(", est[2], "-", est[3], ")^a")
  } 
  
  #pub.bias <- ies.da %>% 
  #  mutate(y = yi/sqrt(vi), x = 1/sqrt(vi)) %>% 
  #  lm(y ~ x, data = .) %>% 
  #  summary()  
  
  #using metaprob from package meta
  # calculate egger's test
  if (nrow(da) <= 2){
    egger.t <- "-"
    egger.pval <- "-"
  } else{
    meta.prop <- metaprop(event = unlist(da[,name]),
      n = da$NumbervalidateNTM,
      method = "Inverse",
      sm = "PFT")
    pub.bias <- metabias(meta.prop, method.bias = "Egger", k.min = 2)
    egger.t <- format(round(pub.bias$statistic,1), nsmall = 1)
    egger.pval <- format_pval(pub.bias$pval)
  }
  strain_res <- rbind(strain_res, 
                      data.frame(NTM_species = name, 
                                 num_study = num_study,
                                 prev_NTM = strain.prev,
                                 n_N = n_N,
                                 I2 = I2,
                                 Q_pval = Q.pval,
                                 egger.t = egger.t,
                                 egger.pval = egger.pval,
                                 est = pes[[1]]
                                 )
                )
}

#write table
write.xlsx(strain_res, file = "Table 3.NTM strain distribution in China.xlsx", asTable = FALSE, overwrite = TRUE)

tab <- strain_res %>% 
  select(NTM_species, num_study, prev_NTM, n_N, I2, Q_pval, egger.t, egger.pval) %>%
  mutate_if(is.character, stri_enc_toutf8)
#for(i in 1:nrow(tab)){
#  tab[i,2] <- r2rtf::utf8Tortf(tab[i,2])
#}	

fileName <- "Table 4.NTM strain distribution in South and North China.rtf"
tabName <- "Table 4.NTM strain distribution in South and North China"
col_header <- "NTM Type|Number of Study|NTM ProportIon|n/N|Cross Study Heterogeneity|Heterogeneity p value|Egger's test t|Egger's test p value"
rel_width <- c(1,1,3,2,1,1,1,1)

tbl <- tab %>%
  # Table title
  rtf_title(
    r2rtf::utf8Tortf(stri_enc_toutf8(tabName))
  ) %>%
  # First row of column header 
  rtf_colheader(r2rtf::utf8Tortf(stri_enc_toutf8(col_header)),
                col_rel_width = rel_width,
                text_background_color = "snow3",
                border_top = rep("single",8)
  ) %>%
  # Table body
  rtf_body(
    col_rel_width = rel_width,
    text_justification = c("l", rep("c", 7)),
    border_top = rep("single",8),
    border_bottom = "",
    border_last="single",
    #group_by = c("TESTNAME"),
    #group_by = "NTM_cat",
    cell_vertical_justification = c("merge_rest",rep("middle",7))
  ) %>%
  rtf_footnote(
    footnote=r2rtf::utf8Tortf(stri_enc_toutf8("{^a}, 95% CI by Clopper-Pearson exact method")),
    as_table = F		   
  )   

#tbls <- list(tbl_1, tbl_2)
tbl %>% 
  rtf_encode() %>%
  # Save to a file
  write_rtf(fileName)  

ies.da <- escalc(xi = ntm$SGM, ni = ntm$NumbervalidateNTM, add = 0, measure = "PFT")
da <- cbind(ntm, ies.da) %>%
  filter(NumbervalidateNTM > 0) %>%
  filter(SGM > 0)

#doing subgroup analysis using metaprop
pes.summary <- metaprop(event = SGM, n = NumbervalidateNTM,
                        data = da, studlab = Study, sm = "PFT",
                        method.tau = "REML", method.ci = "WS",
                        backtransf = T, pscale = 1, 
                        subgroup = RegionSN, subgroup.name = "Region",
                        print.subgroup.name  = TRUE, sep.subgroup = ":",
                        common = FALSE, random = TRUE, tau.common = F,
                        test.subgroup = TRUE, overall.hetstat = TRUE, 
                        overall = TRUE) 
# use weight similar to rma, namely ni = 1/se^2
pes.summary$n.harmonic.mean <- 1/pes.summary$seTE^2
pes.summary$n.harmonic.mean.w <- 1/pes.summary$seTE.random.w^2
pes.summary$n.harmonic.mean.ma <- 1/pes.summary$seTE.random^2

png(file = "Figure 4. SGM isolation rate.png", width = 12, height = 21,units = "in", res = 600)
forest(pes.summary, 
       sortvar = TE,
       common = FALSE, 
       print.tau2 = TRUE, 
       print.Q = TRUE, 
       print.pval.Q = TRUE, 
       print.I2 = TRUE, 
       rightcols = FALSE, 
       pooled.totals = FALSE, 
       weight.study = "random", 
       hetstat = TRUE,
       leftcols = c("studlab", "event", "n", "effect", "ci"), 
       leftlabs = c("Study", "SGM isolate", "Total", "Isolation rate (%)", "95% C.I."), 
       xlab = "Isolation rate of SGM (%)", 
       smlab = "", 
       xlim = c(0,100), 
       pscale = 100, 
       squaresize = 0.5, 
       fs.hetstat = 10, 
       digits = 2, 
       col.square = "navy", 
       col.square.lines = "navy", 
       col.diamond = "maroon", 
       col.diamond.lines = "maroon")
dev.off()  

ies.da <- escalc(xi = ntm$RGM, ni = ntm$NumbervalidateNTM, add = 0, measure = "PFT")
da <- cbind(ntm, ies.da) %>%
  filter(NumbervalidateNTM > 0) %>%
  filter(RGM > 0)
#doing subgroup analysis using metaprop
pes.summary <- metaprop(event = RGM, n = NumbervalidateNTM,
                        data = da, studlab = Study, sm = "PFT",
                        method.tau = "REML", method.ci = "WS",
                        backtransf = T, pscale = 1, 
                        subgroup = RegionSN, subgroup.name = "Region",
                        print.subgroup.name  = TRUE, sep.subgroup = ":",
                        common = FALSE, random = TRUE, tau.common = F,
                        test.subgroup = TRUE, overall.hetstat = TRUE, 
                        overall = TRUE) 
# use weight similar to rma, namely ni = 1/se^2
pes.summary$n.harmonic.mean <- 1/pes.summary$seTE^2
pes.summary$n.harmonic.mean.w <- 1/pes.summary$seTE.random.w^2
pes.summary$n.harmonic.mean.ma <- 1/pes.summary$seTE.random^2
png(file = "Figure 5. RGM isolation.png", width = 12, height = 21,units = "in", res = 600)
forest(pes.summary, 
       sortvar = TE,
       common = FALSE, 
       print.tau2 = TRUE, 
       print.Q = TRUE, 
       print.pval.Q = TRUE, 
       print.I2 = TRUE, 
       rightcols = FALSE, 
       pooled.totals = FALSE, 
       weight.study = "random", 
       hetstat = TRUE,
       leftcols = c("studlab", "event", "n", "effect", "ci"), 
       leftlabs = c("Study", "RGM isolate", "Total", "Isolation rate (%)", "95% C.I."), 
       xlab = "Isolation rate of RGM (%)", 
       smlab = "", 
       xlim = c(0,100), 
       pscale = 100, 
       squaresize = 0.5, 
       fs.hetstat = 10, 
       digits = 2, 
       col.square = "navy", 
       col.square.lines = "navy", 
       col.diamond = "maroon", 
       col.diamond.lines = "maroon")
dev.off()  

