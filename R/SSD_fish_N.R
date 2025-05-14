library(readxl)
library(rcompanion)
# library(tidyr)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(stringr)
library(hydroGOF)
# library(GoFKernel)

setwd("C:/Users/DELL/Desktop/zhouj8/01Research/04Global_CF/04data/fish_spicies/occ2range4fish-master/")
# 

# 
# # output
out_fig <- "map/"


# data processing ====
# function-calculate PDF (potentially-disappeared fraction of species) for fish
PDF_fish <- function(pair_survival){
  # maximal survival species in a ecoregion/biome
  max_survival <- pair_survival %>%
    group_by(ID) %>%
    mutate(max = max(Species_richness)) %>%
    ungroup()  
  # calculate PDF
  pair_survival <- max_survival %>%
    mutate(PDFs = 1-Species_richness/max) %>%
    select(-max)

  return(pair_survival)}

# logistic regression ====
# nls.control(maxiter = 1000000000)

# ssd_fun <- function(x,a,b) {1 / (1 + a*exp (b*x))}

# ssd_fun_inv <- function(x,a,b) {1 / (1 + exp ((a-x)/b))} # assume a=c50
ssd_fun <- function(x,a,b) {1 / (1 + exp ((a-x)/b))} # assume a=c50
ssd_log_fun <- function(x,a,b) {1 / (1 + exp ((a-log10(x))/b))} # assume a=c50
# Pdfs =1 / (1 + exp (b*(a-x)))

# estimate c50
# nls regress PDFs and N concentration
findc_fun <- function(x,y){
  loc.close <- which(abs(y-0.5)==min(abs(y-0.5)))
  ya <- y[loc.close]
  # only one dot of closing to c50
  if(length(ya)==1){
    if(ya-0.5==0){x50 <- x[loc.close]}
    else{
      xa <- x[loc.close]
      if(ya-0.5 > 0 ){
        if(identical(y[loc.close+1], numeric(0))){
          print("numeric(y[loc.close+1])")
          loc.close <- loc.close-1
          yb <- y[loc.close]
          xb<- x[loc.close]}
        else{
          loc.close <- loc.close+1
          yb <- y[loc.close]
          xb<- x[loc.close]}
        }
      else{
        if(identical(y[loc.close-1], numeric(0))){
          loc.close <- loc.close+1
          yb <- y[loc.close]
          xb<- x[loc.close]}
        else{loc.close <- loc.close-1
        yb <- y[loc.close]
        xb<- x[loc.close]}
      }
      x50 <- (xa-xb)/(ya-yb)*(0.5-ya)+xa}
  }
  # getting two dots of closing c50
  else{
    ya1 <- ya[1]
    ya2 <- ya[2]
    xa <- x[loc.close]
    xa1 <- xa[1]
    xa2 <- xa[2]
    x50 <- (xa1-xa2)/(ya1-ya2)*(0.5-ya1)+xa1
  }
  return(log10(x50))
}

# regression
SSD_fit_fun <- function(Nconc,PDF,weight,regionid,n,nlc=nls.control(maxiter = 100000)){
  # regress with observed c50
  a.0 <- findc_fun(Nconc,PDFs)
  # start <- list(a=a.0, b=1.5086537) # for Nconc
  start <- list(a=a.0, b=0.1) # for log10 Nconc
  t1 <- try(nls(PDFs ~ ssd_fun(log10(Nconc), a, b),weights = weight, control = nlc, start = start))

  # regress with linear-estimated a and b
  SSD_fit.0 <- lm(log(1/(PDFs+0.0000001)-0.99999999) ~ Nconc)  # for Nconc

  start <- list(a=-coef(SSD_fit.0)[1]/coef(SSD_fit.0)[2], b=-1/coef(SSD_fit.0)[2])
  t2 <- try(nls(PDFs ~ ssd_fun(log10(Nconc), a, b),weights = weight,control = nlc, start = start))
  pseudo_r2_t2<- pseudo_r2_fun(Nconc, PDFs, t2)
  if(inherits(t1, "try-error")) {
    print(paste("Ecoregion add", n, "weight:", regionid, "has", length(Nconc),
                "PDFs-Nconc data,but encounter singular gradient of c50",sep=" "))
    SSD_fit <-t2
    if(inherits(t2, "try-error")) {
      # print(paste("Ecoregion add", n, "weight:", regionid, "has", length(Nconc),
      #             "PDFs-Nconc data,but encounter singular gradient of lm",sep=" "))
    }
  }
  else {
    pseudo_r2_t1<- pseudo_r2_fun(Nconc, PDFs, t1)
    if(inherits(t2, "try-error")){
      SSD_fit <- t1
    }
    else{
      if(pseudo_r2_t1<pseudo_r2_t2){
        SSD_fit <- t2
      }
      else{SSD_fit <- t1}}
  }
  SSD_fit
}


# pseudo-R2 ====

nullfunct <- function(x, m) {m}
# pseudo-R2
pseudo_r2_fun <- function(x, y, fit) {
  null_model <- nls(y ~ nullfunct(x, m), start = list(m = 0.5))
  pseudo_r2 <- nagelkerke(fit, null = null_model)
  pseudo_r2 <- pseudo_r2$Pseudo.R.squared.for.model.vs.null["Cox and Snell (ML)",]
  pseudo_r2
}
# McFadden's pseudo-R2
# pseudo_r2_McFadden_fun <- function(x, y, fit) {
#   null_model <- nls(y ~ nullfunct(x, m), start = list(m = 1))
#   pseudo_r2_McFadden <- 1-logLik(fit)/logLik(null_model)
#   pseudo_r2_McFadden
# }

# function- select PDF and N concentration
NRMSE_fun<- function(sim,obs){
  error = sim-obs
  squaredError = error*error
  NRMSE = sqrt(sum(squaredError) / length(squaredError))/mean(obs)
  NRMSE
}

# Plot function and pseudo r2
lm_eqn <- function(a,b,r2){
  # eq <- substitute(italic(y) == frac(1,1+exp(frac(a-log[10]~~italic(x), b)))* ", Pseudo"~~italic(r)^2~"="~r2, 
  #                  list(a = format(unname(a), digits = 2),
  #                       b = format(unname(b), digits = 2),
  #                       r2 = format(unname(r2), digits = 2)))
  eq <- substitute("Pseudo"~~r^2~"="~r2, 
                   list(r2 = format(unname(r2), digits = 2)))
  as.character(as.expression(eq))
}

#### calculate PDF for ecoregions and biomes

pair_ecoregion_survival_N <- PDF_fish(pair_ecoregion_survival)
pair_biome_survival_N <- PDF_fish(pair_biome_survival)
pair_realm_survival_N <- PDF_fish(pair_realm_survival)
pair_habitat_survival_N <- PDF_fish(pair_habitat_survival)
#### calibration function
#### Calibration ecoregion and biome
#### create criteria table before run the loop
criteria_biome_N <- as.data.frame(unique(pair_biome_survival_N$ID)) %>%
  mutate(Pseudo_r2=NaN, NRMSE_mean_PDF=NaN, NRMSE_std_PDF=NaN,PDF0= NaN,a=NaN, b=NaN)
names(criteria_biome_N)[1] <-'ID'

criteria_ecoregion_N <- as.data.frame(unique(pair_ecoregion_survival_N$ID)) %>%
  mutate(Pseudo_r2=NaN, NRMSE_mean_PDF=NaN, NRMSE_std_PDF=NaN,PDF0= NaN,a=NaN, b=NaN)
names(criteria_ecoregion_N)[1] <-'ID'

criteria_realm_N <- as.data.frame(unique(pair_realm_survival_N$ID)) %>%
  mutate(Pseudo_r2=NaN, NRMSE_mean_PDF=NaN, NRMSE_std_PDF=NaN,PDF0= NaN,a=NaN, b=NaN)
names(criteria_realm_N)[1] <-'ID'

criteria_habitat_N <- as.data.frame(unique(pair_habitat_survival_N$ID)) %>%
  mutate(Pseudo_r2=NaN, NRMSE_mean_PDF=NaN, NRMSE_std_PDF=NaN,PDF0= NaN,a=NaN, b=NaN)
names(criteria_habitat_N)[1] <-'ID'


# call loop for plot and criteria # shift+ctrl+F10 to restart R if plot cannot come out
# for biome:
pair_survival <- pair_biome_survival_N
criteria <- criteria_biome_N
# for ecoregion: 
pair_survival <- pair_ecoregion_survival_N
criteria <- criteria_ecoregion_N
# for realm:
pair_survival <- pair_realm_survival_N
criteria <- criteria_realm_N
# for Major habitat: 
pair_survival <- pair_habitat_survival_N
criteria <- criteria_habitat_N
# calibration_fun <- function(pair_survival){}

for(regionid in unique(pair_survival$ID)){
  loc <- which(pair_survival[,'ID']== regionid)
  
  # model calibration
  pair_survival_temp <-pair_survival[loc,]
  
  # skip those ecoregion with too few data
  
  if(length(pair_survival_temp$Nconc)<=3){
    print(paste("Ecoregion:", regionid, "has", length(pair_survival_temp$Nconc),
                "PDFs-Nconc data.",sep=" "))
    next
  }
  # select N concentration
  Nconc <- dplyr::pull(pair_survival_temp,Nconc) %>%
    append(.,0)
  # select PDF of a region/biome
  PDFs <- dplyr::pull(pair_survival_temp,PDFs) %>%
    append(.,0)
  ##### estimate a,b for nls regression
  # add weighting coefficients
  weight = rep(1, length(Nconc))
  # estimate c50 and regress with b(c50-x)
  SSD_fit <- SSD_fit_fun(Nconc,PDF,weight,regionid,n="no",nlc=nls.control(maxiter = 10000000))
  # estimate c50 and regress with b(c50-x)
  
#   # add weights if PDF0 >0.01
  PDF0 <- try(predict(SSD_fit,data.frame(Nconc=c(0))))
  if(inherits(PDF0, "try-error")){
    print(paste("Ecoregion:", regionid, "skip",sep=" "))
    next
  }
  else{
    # print(paste("Ecoregion:", regionid, "fit",sep=" "))
  #   if(PDF0 <=0.01){
  #   SSD_fit <-SSD_fit
  # }
  # else{
  #   weight[length(Nconc)] = 10*length(Nconc)
  #   SSD_fit <- SSD_fit_fun(Nconc,PDF,weight,regionid,n="10L",nlc=nls.control(maxiter = 100000))
  #   PDF0 <- predict(SSD_fit,data.frame(Nconc=c(0)))
  #   # whether PDF0 <0.01 after adding 10L weight
  #   if(PDF0 >0.01){
  #     weight[length(Nconc)] = 100*length(Nconc)
  #     SSD_fit <- SSD_fit_fun(Nconc,PDF,weight,regionid,n="100L",nlc=nls.control(maxiter = 100000))
  #     PDF0 <- predict(SSD_fit,data.frame(Nconc=c(0)))
  #     if(PDF0 >0.01){
  #       weight[length(Nconc)] = 500*length(Nconc)
  #       SSD_fit <- SSD_fit_fun(Nconc,PDF,weight,regionid,n="500L",nlc=nls.control(maxiter = 100000))
  #       PDF0 <- predict(SSD_fit,data.frame(Nconc=c(0)))}
  #   }
  # }
  
  
  # remove Pseudo point(0,0)
  Nconc <-Nconc[1:length(Nconc)]
  PDFs <-PDFs[1:length(PDFs)]
  
  sim_PDF <- predict(SSD_fit,data.frame(Nconc=Nconc))
  # pseudo-R2 and NRMSE ====
  pseudo_r2 <- pseudo_r2_fun(Nconc, PDFs, SSD_fit)
  NRMSE_mean <- NRMSE_fun(sim_PDF,PDFs) 
  NRMSE_std <- nrmse(sim_PDF,PDFs)/100

  
  parameter<- c(unlist(summary(SSD_fit)['parameters']))
  criteria[which(criteria[,'ID']== regionid),'Pseudo_r2'] <- pseudo_r2
  criteria[which(criteria[,'ID']== regionid),'NRMSE_mean_PDF'] <- NRMSE_mean
  criteria[which(criteria[,'ID']== regionid),'NRMSE_std_PDF'] <- NRMSE_std
  criteria[which(criteria[,'ID']== regionid),'PDF0'] <- PDF0
  criteria[which(criteria[,'ID']== regionid),'a'] <- parameter[1]
  criteria[which(criteria[,'ID']== regionid),'b'] <- parameter[2]
  
  
  # #### log nls
  # lgNconc <- log(Nconc)
  # lgPDFs <- log(PDFs)
  # 
  # # estimate c50 and regress with c50+bx
  # a.0 <- findc_fun(lgNconc,lgPDFs)
  # start <- list(a=a.0, b=1.5086537)
  # t <- try(SSD_fit <- nls(lgPDFs ~ ssd_fun(lgNconc, a, b),weights = weight, control = nlc, start = start))
  # # run with estimating a+bx if c50+ bx calls error of singular gradient
  # if(inherits(t, "try-error")) {
  #   print(paste("Ecoregion:", regionid, "has", length(pair_survival_temp$Nconc),
  #               "PDFs-Nconc data,but encounter singular gradient of c50",sep=" "))
  #   # lg(1/pdfs-1) = b(a-x)
  #   
  #   SSD_fit.0 <- lm(log(1/(lgPDFs+0.0001)-0.99) ~ lgNconc)
  #   start <- list(a=-coef(SSD_fit.0)[1]/coef(SSD_fit.0)[2], b=-coef(SSD_fit.0)[2])
  #   SSD_fit <- nls(lgPDFs ~ ssd_fun(lgNconc, a, b),weights = weight,control = nlc, start = start)
  # }
  # else {
  #   SSD_fit <- t
  # }
  # 
  # 
  # # pseudo-R2 ====
  # Pseudo_r2_log2 <-pseudo_r2_fun(lgNconc, lgPDFs, SSD_fit)
  # parameter<- c(unlist(summary(SSD_fit)['parameters']))
  # criteria[which(criteria[,'ID']== regionid),'Pseudo_r2_log2'] <- pseudo_r2_log2
  # criteria[which(criteria[,'ID']== regionid),'loga'] <- parameter[1]
  # criteria[which(criteria[,'ID']== regionid),'logb'] <- parameter[2]

  # plot regression
  df_region <- as.data.frame(cbind(Nconc,PDFs))
  plot1<-ggplot(df_region,aes(Nconc, PDFs))+
    geom_point(aes(colour='Observed'))+ # c(0.8, 0.7)
    stat_function(fun = ssd_log_fun, args = list(coef(SSD_fit)["a"], coef(SSD_fit)["b"]),aes(colour='Predicted'))+
    scale_color_manual(values = c('#DC582A','#006400'),
                       labels = c("Observed",'Predicted'),
                       guide = guide_legend(override.aes = list(
                         shape = c(19, NA),
                         linetype = c("blank", "solid")))) + 
    scale_shape_manual(values = 18:21) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.1)) +
    labs(title = str_wrap(regionid, 35), x=expression('N concentration (mg/L)'), y='PDF (-)',
         color = "") +
    theme(plot.title = element_text(size = rel(2),hjust = 0.5),
          axis.title = element_text(size = rel(2)),
          axis.text = element_text(size = rel(2)),
          axis.line = element_line(color="grey50"),
          legend.text = element_text(size = rel(2)),
          legend.title =element_text(size = rel(2)),
          legend.position=c(0.65, 0.3),
          plot.margin = margin(t = 5, r = 10, b = 5, l = 5, unit = "pt"),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          legend.key = element_rect(fill = NA))+
    annotate(geom="text",
             label = lm_eqn(unlist(parameter[1]),
                             unlist(parameter[2]),
                             unlist(pseudo_r2)),x=0.05,y=1,hjust = 0, vjust = 0,
             size = 8,
             parse = TRUE)
    
  # path and name for saving figures
  out_fig.sigmoid <- paste(out_fig, str_trim(regionid), ".tif", sep='')
  # Save the plot
  tiff(out_fig.sigmoid,width=900, height=1000,res=150,compression = "lzw")
  print(plot1)
  dev.off()}
  
  
}
# for(regionid in unique(pair_survival$ID)){
#   loc <- which(pair_survival[,'ID']== regionid)
#   
#   # model calibration
#   pair_survival_temp <-pair_survival[loc,]
#   
#   # skip those ecoregion with too few data
#   
#   if(length(pair_survival_temp$Nconc)<=3){
#     print(paste("Ecoregion:", regionid, "has", length(pair_survival_temp$Nconc),
#                 "PDFs-Nconc data.",sep=" "))
#     next
#   }
#   # select N concentration
#   Nconc <- dplyr::pull(pair_survival_temp,Nconc) %>%
#     append(.,0)
#   # select PDF of a region/biome
#   PDFs <- dplyr::pull(pair_survival_temp,PDFs) %>%
#     append(.,0)
#   ##### estimate a,b for nls regression
#   # add weighting coefficients
#   weight = rep(1, length(Nconc))
#   # estimate c50 and regress with b(c50-x)
#   SSD_fit <- SSD_fit_fun(Nconc,PDF,weight,regionid,n="no",nlc=nls.control(maxiter = 10000000))
#   # estimate c50 and regress with b(c50-x)
#   
#   #   # add weights if PDF0 >0.01
#   PDF0 <- try(predict(SSD_fit,data.frame(Nconc=c(0))))
#   if(inherits(PDF0, "try-error")){
#     print(paste("Ecoregion:", regionid, "skip",sep=" "))
#     next
#   }
#   else{
#     # print(paste("Ecoregion:", regionid, "fit",sep=" "))
#     #   if(PDF0 <=0.01){
#     #   SSD_fit <-SSD_fit
#     # }
#     # else{
#     #   weight[length(Nconc)] = 10*length(Nconc)
#     #   SSD_fit <- SSD_fit_fun(Nconc,PDF,weight,regionid,n="10L",nlc=nls.control(maxiter = 100000))
#     #   PDF0 <- predict(SSD_fit,data.frame(Nconc=c(0)))
#     #   # whether PDF0 <0.01 after adding 10L weight
#     #   if(PDF0 >0.01){
#     #     weight[length(Nconc)] = 100*length(Nconc)
#     #     SSD_fit <- SSD_fit_fun(Nconc,PDF,weight,regionid,n="100L",nlc=nls.control(maxiter = 100000))
#     #     PDF0 <- predict(SSD_fit,data.frame(Nconc=c(0)))
#     #     if(PDF0 >0.01){
#     #       weight[length(Nconc)] = 500*length(Nconc)
#     #       SSD_fit <- SSD_fit_fun(Nconc,PDF,weight,regionid,n="500L",nlc=nls.control(maxiter = 100000))
#     #       PDF0 <- predict(SSD_fit,data.frame(Nconc=c(0)))}
#     #   }
#     # }
#     
#     
#     # remove Pseudo point(0,0)
#     Nconc <-Nconc[1:length(Nconc)]
#     PDFs <-PDFs[1:length(PDFs)]
#     
#     sim_PDF <- predict(SSD_fit,data.frame(Nconc=Nconc))
#     # pseudo-R2 and NRMSE ====
#     pseudo_r2 <- pseudo_r2_fun(Nconc, PDFs, SSD_fit)
#     NRMSE_mean <- NRMSE_fun(sim_PDF,PDFs) 
#     NRMSE_std <- nrmse(sim_PDF,PDFs)/100
#     
#     
#     parameter<- c(unlist(summary(SSD_fit)['parameters']))
#     criteria[which(criteria[,'ID']== regionid),'Pseudo_r2'] <- pseudo_r2
#     criteria[which(criteria[,'ID']== regionid),'NRMSE_mean_PDF'] <- NRMSE_mean
#     criteria[which(criteria[,'ID']== regionid),'NRMSE_std_PDF'] <- NRMSE_std
#     criteria[which(criteria[,'ID']== regionid),'PDF0'] <- PDF0
#     criteria[which(criteria[,'ID']== regionid),'a'] <- parameter[1]
#     criteria[which(criteria[,'ID']== regionid),'b'] <- parameter[2]
#     
#     
#     # #### log nls
#     # lgNconc <- log(Nconc)
#     # lgPDFs <- log(PDFs)
#     # 
#     # # estimate c50 and regress with c50+bx
#     # a.0 <- findc_fun(lgNconc,lgPDFs)
#     # start <- list(a=a.0, b=1.5086537)
#     # t <- try(SSD_fit <- nls(lgPDFs ~ ssd_fun(lgNconc, a, b),weights = weight, control = nlc, start = start))
#     # # run with estimating a+bx if c50+ bx calls error of singular gradient
#     # if(inherits(t, "try-error")) {
#     #   print(paste("Ecoregion:", regionid, "has", length(pair_survival_temp$Nconc),
#     #               "PDFs-Nconc data,but encounter singular gradient of c50",sep=" "))
#     #   # lg(1/pdfs-1) = b(a-x)
#     #   
#     #   SSD_fit.0 <- lm(log(1/(lgPDFs+0.0001)-0.99) ~ lgNconc)
#     #   start <- list(a=-coef(SSD_fit.0)[1]/coef(SSD_fit.0)[2], b=-coef(SSD_fit.0)[2])
#     #   SSD_fit <- nls(lgPDFs ~ ssd_fun(lgNconc, a, b),weights = weight,control = nlc, start = start)
#     # }
#     # else {
#     #   SSD_fit <- t
#     # }
#     # 
#     # 
#     # # pseudo-R2 ====
#     # Pseudo_r2_log2 <-pseudo_r2_fun(lgNconc, lgPDFs, SSD_fit)
#     # parameter<- c(unlist(summary(SSD_fit)['parameters']))
#     # criteria[which(criteria[,'ID']== regionid),'Pseudo_r2_log2'] <- pseudo_r2_log2
#     # criteria[which(criteria[,'ID']== regionid),'loga'] <- parameter[1]
#     # criteria[which(criteria[,'ID']== regionid),'logb'] <- parameter[2]
#     
#     # plot regression
#     df_region <- as.data.frame(cbind(Nconc,PDFs))
#     plot1<-ggplot(df_region,aes(Nconc, PDFs))+
#       geom_point(aes(colour='Observed'))+ # c(0.8, 0.7)
#       stat_function(fun = ssd_log_fun, args = list(coef(SSD_fit)["a"], coef(SSD_fit)["b"]),aes(colour='Predicted'))+
#       scale_color_manual(values = c('#DC582A','#006400'),
#                          labels = c("Observed",'Predicted'))+
#       scale_shape_manual(values = 18:21) + 
#       scale_x_continuous(expand = c(0, 0)) +
#       scale_y_continuous(expand = c(0, 0), limits = c(0,1.2)) +
#       labs(title = str_wrap(regionid, 20), x=expression('N concentration (mg/L)'), y='PDF (-)',colour = "Legend")+
#       theme(plot.title = element_text(size = rel(3),hjust = 0.5),
#             axis.title = element_text(size = rel(3)),
#             axis.text = element_text(size = rel(3)),
#             axis.line = element_line(color="grey50"),
#             legend.text = element_text(size = rel(2)),
#             legend.title =element_text(size = rel(2)),
#             legend.position=c(0.8, 0.15),
#             plot.margin = margin(t = 5, r = 10, b = 5, l = 5, unit = "pt"),
#             panel.background = element_rect(fill = "white", colour = "grey50"))+
#       annotate(geom="text",
#                label = lm_eqn(unlist(parameter[1]),
#                               unlist(parameter[2]),
#                               unlist(pseudo_r2)),x=0.05,y=1,hjust = 0, vjust = 0,
#                size = 12,
#                parse = TRUE)
#     
#     # path and name for saving figures
#     out_fig.sigmoid <- paste(out_fig,regionid,".tif",sep='')
#     # Save the plot
#     tiff(out_fig.sigmoid,width=900, height=1000,res=150,compression = "lzw")
#     print(plot1)
#     dev.off()}
#   
#   
# }
# criteria
criteria_biome_N <- criteria
write.csv(criteria_biome_N,file = paste(out_fig,'/criteria_biome_(c50-logx)b.csv',sep=''),na="NA")

criteria_ecoregion_N <- criteria
write.csv(criteria_ecoregion_N,file = paste(out_fig,'/criteria_ecoregion_(c50-logx)b.csv',sep=''),na="NA")

criteria_realm_N <- criteria
write.csv(criteria_realm_N,file = paste(out_fig,'/criteria_realm_b(c50-x).csv',sep=''),na="NA")

criteria_habitat_N <- criteria
write.csv(criteria_habitat_N,file = paste(out_fig,'/criteria_habitat_b(c50-x).csv',sep=''),na="NA")

# grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
dev.off()
dev.new()
  
