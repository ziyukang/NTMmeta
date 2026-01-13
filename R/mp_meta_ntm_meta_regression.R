# meta analysis of prevalence of NTM species 

library(tidyverse)
library(metafor) 
library(readxl)
library(r2rtf)
library(stringi) # 字符串处理
require(openxlsx)
require(meta)
require(DescTools)

format_pval <- function(pval){
  
  if(pval < 0.001) pval <- "< 0.001"
  else pval <- format(round(pval,3), nsmall = 3)
  return(pval)
}

ntm <- read_excel("MP meta分析数据收集表核查prev_eng.xlsx",skip = 1) %>%
  select(FirstAuthor0, FirstAuthor, PublicationYear, StudyDesign, 
         StudyRegion, Regions, ScreeningPatientNumber, 
         prevalence, 总计) %>%
  mutate(NumbervalidateNTM = as.numeric(总计)) %>%
  mutate(Regions = factor(Regions, levels = c("东北","西北","中部","东南","西南"))) %>%
  filter(NumbervalidateNTM > 0) %>%
  filter(prevalence == "是" & !is.na(NumbervalidateNTM)) %>%
  mutate(Study = paste(FirstAuthor, PublicationYear, sep = ",")) %>%
  mutate(ScreeningPatientNumber = as.numeric(ScreeningPatientNumber))

 
  #纳入所有21个研究，排除刘超(2023)并绘制森林图
  #make forest plot for NTM患病率
  da <- ntm %>%
    filter(!FirstAuthor0 == "刘超")
  #da$Study <- paste(da$FirstAuthor, da$PublicationYear, sep = ",")
  #Only choose one of the three transformation methods 
  ies.da <- escalc(xi = NumbervalidateNTM, ni = ScreeningPatientNumber, data = da, add = 0, digits = 10, measure = "PFT")
  ies.da <- cbind(da, ies.da)

  m_reg <- rma(yi, vi, data =ies.da, mods =~ PublicationYear, method = "REML", weighted = TRUE)
  m_reg
  ### draw plot
  mytransf <- function(x, targs) {
    round(transf.ipft.hm(x, targs) * 100,0)
  }
  transf=transf.ipft.hm, targ = list(ni=1/(pes.da$se)^2), 
  png(file = "Figure S3. NTM isolation rate meta regression against publication years(21个研究，排除刘超(2023)).png", width = 9, height = 10,units = "in", res = 300)
  regplot(m_reg, xlim=c(2016,2025), ylim = c(0,0.5), xlab="Publication Year", 
          ylab = "NTM detection rate",
          digits=1, las=1, bty="l",
          offset=c(1.6,0.8), labsize=0.9, psize = 0.7, pch = 2)
  text(2016, 0.3, pos = 4, paste0("y = ", round(m_reg$b[1],4), "+", round(m_reg$b[2],4), "*x"))
  text(2016, 0.27, pos = 4, paste0("p value for slope: ", round(m_reg$pval[2],4)))
  dev.off()
  pes.summary <- metaprop(event = NumbervalidateNTM, n = ScreeningPatientNumber,
                          data = ies.da, studlab = Study, sm = "PFT",
                          method.tau = "REML", method.ci = "WS",
                          backtransf = TRUE, pscale = 1,
                          common = FALSE, random = TRUE,
                          test.subgroup = TRUE, overall.hetstat = TRUE, 
                          overall = TRUE) 
  png(file = "Figure S2. NTM isolation rate across years(21个研究，排除刘超(2023)).png", width = 9, height = 10,units = "in", res = 300)
  forest(pes.summary, 
         sortvar = PublicationYear,
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
         leftlabs = c("Study", "NTM isolate", "Total", "Isolation rate (%)", "95% C.I."), 
         xlab = "Isolation rate of NTMs (%)", 
         smlab = "", 
         xlim = c(0,50), 
         pscale = 100, 
         squaresize = 0.5, 
         fs.hetstat = 10, 
         digits = 2, 
         col.square = "navy", 
         col.square.lines = "navy", 
         col.diamond = "maroon", 
         col.diamond.lines = "maroon")
  dev.off()  
  