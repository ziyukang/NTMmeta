# meta analysis of prevalence of NTM species 

library(tidyverse)
library(metafor) 
library(readxl)
library(r2rtf)
library(stringi) # 字符串处理
require(openxlsx)
require(meta)
require(DescTools)

set.seed(20251119)
format_pval <- function(pval){
  
  if(pval < 0.001) pval <- "< 0.001"
  else pval <- format(round(pval,3), nsmall = 3)
  return(pval)
}

ntm <- read_excel("MP meta分析数据收集表核查prev_eng.xlsx",skip = 1) %>%
  select(FirstAuthor, PublicationYear, StudyDesign, 
         StudyRegion, Regions, ScreeningPatientNumber, 
         prevalence, NTMpatient, 总计) %>%
  filter(prevalence == "是" & !is.na(NTMpatient)) %>%
  mutate(Study = paste(FirstAuthor, PublicationYear, sep = ",")) %>%
  mutate(ScreeningPatientNumber = as.numeric(ScreeningPatientNumber)) %>%
  mutate(NTMpatient = as.numeric(NTMpatient)) %>%
  mutate(Regions = factor(Regions, levels = c("东南","中部","西南"),
                          labels = c("Southeast","Central","Southwest")))
  
  #PFT transformation methods 
  ies.da <- escalc(xi = unlist(ntm$NTMpatient), ni = ntm$ScreeningPatientNumber, add = 0, digits = 10, measure = "PFT")
  ies.da <- cbind(ntm, ies.da)
  
  ntm_res <- data.frame()
  for(region in c("所有地区",unique(as.character(ntm$Regions)))){
    # estimate across all studies
    if(region == "所有地区"){
      da <- ies.da
      pes.da <- rma(yi, vi, data =da , method = "REML", weighted = TRUE)
      pes <- predict(pes.da, transf = transf.ipft.hm, targ = list(ni=1/(pes.da$se)^2))          
    } else{
      da <- subset(ies.da, Regions == region)
      pes.da <- rma(yi, vi, data =da , method = "REML")
      pes <- predict(pes.da, transf = transf.ipft.hm, targ = list(ni=1/(pes.da$se)^2))          
    }
    
    num_study <- pes.da$k
    n <- sum(da[,"NTMpatient"])
    N <- sum(da$ScreeningPatientNumber)
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
      meta.prop <- metaprop(event = unlist(da[,"NTMpatient"]),
                            n = da$ScreeningPatientNumber,
                            method = "Inverse",
                            sm = "PFT")
      pub.bias <- metabias(meta.prop, method.bias = "Egger", k.min = 2)
      egger.t <- format(round(pub.bias$statistic,1), nsmall = 1)
      egger.pval <- format_pval(pub.bias$pval)
    }
    ntm_res <- rbind(ntm_res, 
                     data.frame(Region = region, 
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
write.xlsx(ntm_res, file = "表1.中国境内非结核分枝杆菌(NTM)病患病率亚组分析.xlsx", asTable = FALSE, overwrite = TRUE)

tab <- ntm_res %>% 
  select(Region, num_study, prev_NTM, n_N, I2, Q_pval, egger.t, egger.pval) %>%
  mutate_if(is.character, stri_enc_toutf8)

for(i in 1: nrow(tab)){
  tab[i,1] = r2rtf::utf8Tortf(stri_enc_toutf8(tab[i,1]))
}

fileName <- "表1.中国境内非结核分枝杆菌(NTM)病患病率亚组分析.rtf"
tabName <- "表1.中国境内非结核分枝杆菌(NTM)病患病率亚组分析"
col_header <- "Region|Number of Study|Prevalence|n/N|Cross Study Heterogeneity|Heterogeneity p value|Egger's test t|Egger's test p"
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
    text_justification = c("c", rep("c", 7)),
    border_top = rep("single",8),
    border_bottom = "",
    border_last="single",
    #group_by = c("TESTNAME"),
    cell_vertical_justification = c("middle",rep("middle",8))
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

#make forest plot for NTM患病率
res <- rma(yi, vi, data =ies.da , method = "REML", weighted = TRUE, slab = Study)

pes.summary <- metaprop(event = NTMpatient, n = ScreeningPatientNumber,
                        data = ies.da, studlab = Study, sm = "PFT",
                        #subset = which(ies.da$Regions == "东南"),
                        method.tau = "REML", method.ci = "WS",
                        backtransf = T, pscale = 1, 
                        subgroup = Regions, subgroup.name = "Region",
                        print.subgroup.name  = TRUE, sep.subgroup = ":",
                        common = FALSE, random = TRUE, tau.common = F,
                        test.subgroup = TRUE, overall.hetstat = TRUE, 
                        overall = TRUE)
# use weight similar to rma, namely ni = 1/se^2
pes.summary$n.harmonic.mean <- 1/pes.summary$seTE^2
pes.summary$n.harmonic.mean.w <- 1/pes.summary$seTE.random.w^2
pes.summary$n.harmonic.mean.ma <- 1/pes.summary$seTE.random^2

# need not to assign func.backtransf, forest.meta use default asim2p function to backtransf
#pes.summary$func.backtransf = "transf.ipft.hm"
#pes.summary$args.backtransf = list(ni = pes.summary$n.harmonic.mean.w)
trace(forest.meta, edit = T)
png(file = "Figure 2. NTM disease prevalence subgroup analysis(9 studies).png", width = 8, height = 7,units = "in", res = 300)
forest(pes.summary, 
       sortvar = TE,
       resid.hetstat = FALSE ,
       common = FALSE, 
       print.tau2 = TRUE, 
       print.Q = TRUE, 
       print.pval.Q = TRUE, 
       print.I2 = TRUE, 
       rightcols = FALSE, 
       pooled.totals = FALSE, 
       weight.study = "random", 
       leftcols = c("studlab", "event", "n", "effect", "ci"), 
       leftlabs = c("Study", "Cases", "Total", "Prevalence(%)", "95% C.I."), 
       xlab = "Prevalence of NTM disease (%)", 
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
untrace(forest.meta)
#identifying outlying studies
ies.da <- escalc(xi = unlist(ntm$NTMpatient), ni = ntm$ScreeningPatientNumber, add = 0, digits = 10, measure = "PFT")
ies.da <- cbind(ntm, ies.da)
da <- ies.da
pes.da <- rma(yi, vi, data =da , slab=Study, method = "REML", weighted = TRUE)

bjplot <- baujat(pes.da, 
                 symbol=19, 
                 xlim=c(0,10), 
                 #ylim = c(0,40),
                 xlab="Contribution to Overall Heterogeneity", 
                 ylab="Influence on Summary Proportion") 
# Label those studies located in the upper quadrants 
bjplot <- bjplot[bjplot$x >= 10 | bjplot$y >= 0.4,] 
text(bjplot$x, bjplot$y, bjplot$slab, pos=1)

inf <- influence(pes.da) 
print(inf, digits=3) 
plot(inf)

# 分地区进行亚组分析，绘制forest plot

#make forest plot for NTM患病率
pes.summary <- metaprop(event = NTMpatient, n = ScreeningPatientNumber,
                        data = da, studlab = Study, sm = "PFT",
                        method.tau = "REML", method.ci = "WS",
                        backtransf = TRUE, pscale = 1,
                        subgroup = Regions, subgroup.name = "Region",
                        print.subgroup.name  = TRUE, sep.subgroup = ":",
                        common = FALSE, random = TRUE,
                        test.subgroup = TRUE, overall.hetstat = TRUE, 
                        overall = TRUE) 
png(file = "NTM病原患病率(9个研究).png", width = 9, height = 10,units = "in", res = 300)
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
       leftlabs = c("Study", "Cases", "Total", "Detection rate", "95% C.I."), 
       xlab = "Detection rate of NTM (%)", 
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

### fit random-effects model
da <- ies.da
res <- rma(yi, vi, data =da , method = "REML", weighted = TRUE, slab = Study)
#res <- predict(pes.da, transf = transf.ipft.hm, targ = list(ni=1/(res$se)^2))

### a little helper function to add Q-test, I^2, and tau^2 estimate info
mlabfun <- function(text, x) {
  list(bquote(paste(.(text)
                    ,
                    " (Q = ", .(fmtx(x$QE, digits=2))
                    ,
                    ", df = ", .(x$k - x$p), ", "
                    ,
                    .
                    (fmtp2(x$QEp)), "; "
                    ,
                    I^2, " = ", .(fmtx(x$I2, digits=1)), "%, "
                    ,
                    tau
                    ^2, " = ", .(fmtx(x$tau2, digits=2)), ")")))
}

mytransf <- function(x, targs) {
  round(transf.ipft.hm(x, targs) * 100,2)
}

### set up forest plot (with 2x2 table counts added; the 'rows' argument is
### used to specify in which rows the outcomes will be plotted)
png(file = "Figure 2. NTM disease prevalence subgroup analysis(9 studies).png", width = 9, height = 10,units = "in", res = 600)
forest(res, xlim=c(-100, 80), at=c(0, 10, 20, 30, 40, 50), 
       xlab = "NTM disease prevalence(%)",# showweights = TRUE,
       transf = mytransf, targs = list(ni=1/(res$se)^2),
       ilab =cbind(NTMpatient, ScreeningPatientNumber), ilab.lab=c("NTM patient","Total patient"),
       ilab.xpos =c(-35,-15), 
       cex=0.75, ylim=c(-3,25), top=4, order= Regions,
       rows =c(3:8,13:14,19), mlab=mlabfun("RE Model for All Studies", res),
       psize = 0.5, header=c("Author(s) and Year","NTM patient prevalence(%)[95% CI]"))

### set font expansion factor (as in forest() above)
op <- par(cex=0.75)

### add additional column headings to the plot
#text(c(-8.75,-5.25), 27, c("Vaccinated", "Control"), font=2)

### add text for the subgroups
text(-100, c(9,15,20), pos=4, c("Southwest"
                                      ,
                                      "Central"
                                      ,
                                      "Southeast"
                                      ), font=4)

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
#res.nw <- rma(yi, vi, subset=(Regions =="西北"), data=da)
#res.nw.pred <- predict(res.nw, transf = mytransf, targ = list(ni=1/(res.nw$se)^2))

#res.ne <- rma(yi, vi, subset=(Regions =="东北"), data=da)
#res.ne.pred <- predict(res.ne, transf = mytransf, targ = list(ni=1/(res.ne$se)^2))

res.ce <- rma(yi, vi, subset=(Regions =="中部"), data=da)
#res.ce.pred <- predict(res.ce, transf = mytransf, targ = list(ni=1/(res.ce$se)^2))

res.sw <- rma(yi, vi, subset=(Regions =="西南"), data=da)
#res.sw.pred <- predict(res.sw, transf = mytransf, targ = list(ni=1/(res.sw$se)^2))

res.se <- rma(yi, vi, subset=(Regions =="东南"), data=da)
#res.se.pred <- predict(res.se, transf = mytransf, targ = list(ni=1/(res.se$se)^2))


### add summary polygons for the three subgroups
addpoly(res.sw, row=17.5, transf = mytransf, targ = list(ni=1/(res.sw$se)^2),mlab=mlabfun("RE Model for Southwest", res.sw))
addpoly(res.ce, row= 11.5, transf = mytransf, targ = list(ni=1/(res.ce$se)^2), mlab=mlabfun("RE Model for Central", res.ce))
addpoly(res.se, row= 1.5, transf = mytransf, targ = list(ni=1/(res.se$se)^2),mlab=mlabfun("RE Model for Southeast", res.se))
#addpoly(res.nw, row= 7.5, transf = mytransf, targ = list(ni=1/(res.nw$se)^2), mlab=mlabfun("RE Model for Northwest", res.nw))
#addpoly(res.se, row= 1.5, transf = mytransf, targ = list(ni=1/(res.ne$se)^2), mlab=mlabfun("RE Model for Northeast", res.ne))

### fit meta-regression model to test for subgroup differences
dat.comp <- data.frame(Regions    = c("东南", "中部", "西南"), 
                       estimate = c(coef(res.se), coef(res.ce), coef(res.sw)), 
                       stderror = c(se(res.se), se(res.ce), se(res.sw)),
                       tau2     = c(res.se$tau2, res.ce$tau2, res.sw$tau2))
res <- rma(estimate, sei = stderror, mods = ~ 0+Regions, method = "FE", data=dat.comp)


### add text for the test of subgroup differences
text(-100, -2, pos=4, cex=0.75, bquote(paste("Test for Subgroup Differences: "
                                             ,
                                             Q
                                             [M], " = ", .(fmtx(res$QM, digits=2))
                                             ,
                                             ", df = ", .(res$p - 1), ", ", .(fmtp2(res$QMp))))) 
dev.off()

