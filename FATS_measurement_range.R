##### This code created by mahdi Abbasi
#### Load requried libraries
rm(list = ls())
library(xlsx)
library(ggplot2)
library(scales)
### Read Excel of evlocity resolution
Ur <- read.xlsx("30 kHz.xlsx",sheetName = "Velocity Resolution")
x_all <- Ur[,1]
U_r <- Ur[,2]
### Modeling start
a_s <- c(seq(0.1,1, length.out =18 ),seq(1,10,length.out = 18),seq(10,100,length.out = 18),seq(100,1000,length.out = 18),
         seq(1000,10000,length.out = 18),seq(10000,100000,length.out = 36))
frequency <- 30
sediment_density <- 2650       ## the density of suspended particle
gamma_coef <- 0.18         ## for quartz particls
speed_sound <- 1500        ## the sound speed in water
omega <- 2*pi*frequency*1000
vescosity <- 1.3e-06      ## the viscosity of water
water_density <- 1000   ## the density of water
zitta_all <- vector("numeric")
zitta_scat <- vector("numeric")
zitta_absor <- vector("numeric")
for (i in 1:126) {
  ##### determine zitta for particles
  ### zitta_scattering form [Sheng and Hay 1988, Richards et al 1996]
  radius_particle <- a_s[i]*10^-06       ## the sphere radius of suspended particles
  acoustic_wavenumber <- (2*pi*frequency*1000)/speed_sound   ## the acoustic wavenumber
  zitta_scattering <- (10*log10(exp(2))/(sediment_density*radius_particle))*((gamma_coef*(acoustic_wavenumber*radius_particle)^4)/(1+(acoustic_wavenumber*radius_particle)^2+(4/3)*(gamma_coef*(acoustic_wavenumber*radius_particle)^4)))
  ### zitta_absorption from Urick 1948
  sigma_sediment <- sediment_density/water_density
  beta <- sqrt(omega/(2*vescosity))
  delta <- 0.5*(1+(9/(2*beta*radius_particle)))
  S <- (9/(4*beta*radius_particle))*(1+(1/beta*radius_particle))
  zitta_absorption <- 10*log10(exp(2))*((acoustic_wavenumber*(sigma_sediment-1)^2)/(2*sediment_density))*(S/(S^2+(sigma_sediment+delta)^2))
  zitta_scat[i] <- zitta_scattering
  zitta_absor[i] <- zitta_absorption
  zitta_all[i] <- (zitta_scattering + zitta_absorption)
 
}



#### ratio_alpaha and SSC
vector_SSC <- c(0.01,0.1,0.2,0.5,1,2,2.75,5.28)
ratio_alpha_SSC <- matrix(NA,nrow = 8, ncol = 126)
for (i in 1:8) {
  for (j in 1:126) {
    ratio_alpha_SSC[i,j] <- zitta_all[j]*vector_SSC[i]
    
  }
  
}

#### plot sound attunation constants as a function of partical size ...
####   Figure 2 of manuscript
plot_zitta <- data.frame(a_s=a_s,zitta_all,zitta_absor,zitta_scat)
tiff("caysi_250.tiff", units="in", width=5, height=5, res=250)
ggplot()+
  geom_line(data=plot_zitta,aes(x=a_s,y=plot_zitta[,2]),col="cyan",size=2)+
  geom_line(data=plot_zitta,aes(x=a_s,y=plot_zitta[,3]),size=1,linetype="solid")+
  geom_line(data=plot_zitta,aes(x=a_s,y=plot_zitta[,4]),size=1.2,linetype="dashed")+
  theme_bw()+
  annotation_logticks()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  xlab("Grain Radius (μm)")+
  ylab("ξ (dB m2 kg-1)")

dev.off()
#### Cross sctional average sound absorbtion coefficient...
####   figure 3 of manuscript
data_plot <- data.frame(a_s=a_s,t(ratio_alpha_SSC))
colnames(data_plot) <- c("a_s","SSC_0.01","SSC_0.1","SSC_0.2","SSC_0.5","SSC_1","SSC_2","SSC_2.39","SSC_5.28")
tiff("SSC_all_250_peak1.tiff", units="in", width=5, height=5, res=250)
ggplot()+
  geom_line(data = data_plot,aes(x=a_s,y=data_plot[,2]),col="sienna4",size=1)+
  geom_line(data = data_plot,aes(x=a_s,y=data_plot[,3]),col="gold1",size=1)+
  geom_line(data = data_plot,aes(x=a_s,y=data_plot[,4]),col="green",size=1)+
  geom_line(data = data_plot,aes(x=a_s,y=data_plot[,5]),col="blue",size=1)+
  geom_line(data = data_plot,aes(x=a_s,y=data_plot[,6]),col="purple",size=1)+
  geom_line(data = data_plot,aes(x=a_s,y=data_plot[,7]),col="#33FFFF",size=1)+
  geom_line(data = data_plot,aes(x=a_s,y=data_plot[,8]),col="red",size=1)+
  geom_line(data = data_plot,aes(x=a_s,y=data_plot[,9]),col="tan4",size=1)+
  theme_bw()+
  annotation_logticks()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(y=expression("<"~alpha~"s>"~"(dB m-1)"),x="Grain Radius (μm)")+
  annotate("text", x = 30000, y = 0.0002, label = "0.01")+
  annotate("text", x = 20000, y = 0.01, label = "0.1")+
  annotate("text", x = 20000, y = 0.02, label = "0.2")+
  annotate("text", x = 20000, y = 0.04, label = "0.5")+
  annotate("text", x = 20000, y = 0.076, label = "1")+
  annotate("text", x = 20000, y = 0.145, label = "2")+
  annotate("text", x = 20000, y = 0.199, label = "2.75")+
  annotate("text", x = 40000, y = 0.38, label = "5.28(kg m-3) ")+
  annotate("text", x = 5, y = 0.22, label = "Viscous absorbtion zone ")+
  annotate("text", x = 1400, y = 0.35, label = "Scattering loss zone")

dev.off()


##### SNR for the peak 1
####  FATS measurement range....
####   Figure 4 of manuscript
SL <- 190
L0 <- 10 
G <- 36
TN <- 120
max_alpha_peak1 <- c(0.000368,0.00368,0.00737,0.01842,0.0368,0.0737,0.1013,0.192)
## SNR_0  SNR for just water without SSC
alpha_0_peak1 <- 0.008
SNR_0_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0_peak1*r-L0+G-TN
}
r_0_peak1 <- seq(100,2331,length.out = 106)
for (i in 1:length(r_0_peak1)) {
  r <- r_0_peak1[i]
  SNR_0_peak1[i] <- SNR(r)
  
}
## SNR_0.01_peak1
alpha_0.01_peak1 <- 0.008+max_alpha_peak1[1]
SNR_0.01_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0.01_peak1*r-L0+G-TN
}
r_0.01_peak1 <- seq(100,2260.5,length.out = 106)
for (i in 1:length(r_0.01_peak1)) {
  r <- r_0.01_peak1[i]
  SNR_0.01_peak1[i] <- SNR(r)
  
}
### SNR_0.1
alpha_0.1_peak1 <- 0.008+max_alpha_peak1[2]
SNR_0.1_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0.1_peak1*r-L0+G-TN
}
r_0.1_peak1 <- seq(100,1501,length.out = 106)
for (i in 1:length(r_0.1_peak1)) {
  r <- r_0.1_peak1[i]
  SNR_0.1_peak1[i] <- SNR(r)
  
}
### SNR_0.2_peak1
alpha_0.2_peak1 <- 0.008+max_alpha_peak1[3]
SNR_0.2_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0.2_peak1*r-L0+G-TN
}
r_0.2_peak1 <- seq(100,1473,length.out = 106)
for (i in 1:length(r_0.2_peak1)) {
  r <- r_0.2_peak1[i]
  SNR_0.2_peak1[i] <- SNR(r)
  
}
### SNR_0.5_peak1
alpha_0.5_peak1 <- 0.008+max_alpha_peak1[4]
SNR_0.5_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0.5_peak1*r-L0+G-TN
}
r_0.5_peak1 <- seq(100,988,length.out = 106)
for (i in 1:length(r_0.5_peak1)) {
  r <- r_0.5_peak1[i]
  SNR_0.5_peak1[i] <- SNR(r)
  
}
#### SNR_1_peak1
alpha_1_peak1 <- 0.008+max_alpha_peak1[5]
SNR_1_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_1_peak1*r-L0+G-TN
}
r_1_peak1 <- seq(100,660.7,length.out = 106)
for (i in 1:length(r_1_peak1)) {
  r <- r_1_peak1[i]
  SNR_1_peak1[i] <- SNR(r)
  
}
#### SNR_2_peak1
alpha_2_peak1 <- 0.008+max_alpha_peak1[6]
SNR_2_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_2_peak1*r-L0+G-TN
}
r_2_peak1 <- seq(100,412,length.out = 106)
for (i in 1:length(r_2_peak1)) {
  r <- r_2_peak1[i]
  SNR_2_peak1[i] <- SNR(r)
  
}
#### SNR_2.75_peak1
alpha_2.75_peak1 <- 0.008+max_alpha_peak1[7]
SNR_2.75_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_2.75_peak1*r-L0+G-TN
}
r_2.75_peak1 <- seq(100,327,length.out = 106)
for (i in 1:length(r_2.75_peak1)) {
  r <- r_2.75_peak1[i]
  SNR_2.75_peak1[i] <- SNR(r)
  
}
### SNR_5.28 peak1
alpha_5.28_peak1 <- 0.008+max_alpha_peak1[8]
SNR_5.28_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_5.28_peak1*r-L0+G-TN
}
r_5.28_peak1 <- seq(100,200,length.out = 106)
for (i in 1:length(r_5.28_peak1)) {
  r <- r_5.28_peak1[i]
  SNR_5.28_peak1[i] <- SNR(r)
  
}
### SNR_12.67 peak1
alpha_12.67_peak1 <- 0.46
SNR_12.67_peak1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_12.67_peak1*r-L0+G-TN
}
r_12.67_peak1 <- 100
for (i in 1:length(r_12.67_peak1)) {
  r <- r_12.67_peak1[i]
  SNR_12.67_peak1[i] <- SNR(r)
  
}

#### SNR_ ggplot for all distance for the first's peak
tiff("SNR_300_peak1.tiff", units="in", width=10, height=5, res=300)
plot(r_0_peak1,SNR_0_peak1,type="l",axes=F,ylim = c(0,50),xlim = c(100,2400),col="black",
     ann=F,xlab = NULL, ylab = NULL,lwd=2)
lines(r_0.01_peak1,SNR_0.01_peak1,type = "l",col="red",lwd=2)
lines(r_0.1_peak1,SNR_0.1_peak1,type = "l",col="blue",lwd=2)
lines(r_0.2_peak1,SNR_0.2_peak1,type = "l",col="green",lwd=2)
lines(r_0.5_peak1,SNR_0.5_peak1,type = "l",col="chocolate",lwd=2)
lines(r_1_peak1,SNR_1_peak1,type = "l",col="deeppink4",lwd=2)
lines(r_2_peak1,SNR_2_peak1,type = "l",col="#999900",lwd=2)
lines(r_2.75_peak1,SNR_2.75_peak1,type = "l",col="darkviolet",lwd=2)
lines(r_5.28_peak1,SNR_5.28_peak1,type = "l",col="#33CCFF",lwd=2)
points(100,SNR_12.67_peak1,col=434,lwd=2)
text(1900,18,"<SSC>=0",srt=-15,cex = 0.8,col = "black",lwd=2)
text(1550,17,"<SSC>=0.01(kg m-1)",srt=-15,cex = 0.8,col = "red",lwd=2)
text(1180,20,"<SSC>=0.1",srt=-20,cex = 0.8,col = "blue",lwd=2)
text(1000,18,"<SSC>=0.2",srt=-26,cex = 0.8,col = "green",lwd=2)
text(820,18,"<SSC>=0.5",srt=-35,cex = 0.8,col = "chocolate",lwd=2)
text(590,18,"<SSC>=1",srt=-50,cex = 0.8,col = "deeppink4",lwd=2)
text(400,18,"<SSC>=2",srt=-65,cex = 0.8,col = "#999900",lwd=2)
text(220,18,"<SSC>=2.75",srt=-75,cex = 0.8,col = "darkviolet",lwd=2)
text(145,18,"<SSC>=5.28",srt=-75,cex = 0.8,col = "#33CCFF",lwd=2)
text(100,18,"<SSC>=12.67",srt=-85,cex = 0.8,col = 434,lwd=2)
axis(side = 4, ylim=c(0,50),col="black",lwd=2,line = -1.5)
mtext(side = 4,text="SNR (dB)",line=0.5)
par(new=T)
plot(x_all,U_r, axes=F, ylim=c(0,max(U_r)), xlab="", ylab="",
     type="l",col="blue", main="",xlim=c(100,2400),lwd=2,lty=2)
axis(side=2,ylim=c(0,max(U_r)),lwd=2,col = "blue",col.axis="blue",line = 1)
mtext(2,text = "Ur (m/s)",line=2.9,col = "blue")
axis(1,pretty(range(100,2400),10))
mtext("Range (m)",side=1,col="black",line=2)
dev.off()

##### SNR for the peak 2
####  FATS measurement range....
####   Figure 4 of manuscript
SL <- 190
L0 <- 10 
G <- 36
TN <- 120
max_alpha <- apply(ratio_alpha_SSC,1,max)
## SNR_0  SNR for just water without SSC
alpha_0 <- 0.008
SNR_0 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0*r-L0+G-TN
}
r_0 <- seq(100,2330,length.out = 106)
for (i in 1:length(r_0)) {
  r <- r_0[i]
  SNR_0[i] <- SNR(r)
  
}
## SNR_0.01
alpha_0.01 <- 0.008+max_alpha[1]
SNR_0.01 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0.01*r-L0+G-TN
}
r_0.01 <- seq(100,2200,length.out = 106)
for (i in 1:length(r_0.01)) {
  r <- r_0.01[i]
  SNR_0.01[i] <- SNR(r)
  
}
### SNR_0.1
alpha_0.1 <- 0.008+max_alpha[2]
SNR_0.1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0.1*r-L0+G-TN
}
r_0.1 <- seq(100,1500,length.out = 106)
for (i in 1:length(r_0.1)) {
  r <- r_0.1[i]
  SNR_0.1[i] <- SNR(r)
  
}
### SNR_0.2
alpha_0.2 <- 0.008+max_alpha[3]
SNR_0.2 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0.2*r-L0+G-TN
}
r_0.2 <- seq(100,1140,length.out = 106)
for (i in 1:length(r_0.2)) {
  r <- r_0.2[i]
  SNR_0.2[i] <- SNR(r)
  
}
### SNR_0.5
alpha_0.5 <- 0.008+max_alpha[4]
SNR_0.5 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_0.5*r-L0+G-TN
}
r_0.5 <- seq(100,683,length.out = 106)
for (i in 1:length(r_0.5)) {
  r <- r_0.5[i]
  SNR_0.5[i] <- SNR(r)
  
}
#### SNR_1
alpha_1 <- 0.008+max_alpha[5]
SNR_1 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_1*r-L0+G-TN
}
r_1 <- seq(100,430,length.out = 106)
for (i in 1:length(r_1)) {
  r <- r_1[i]
  SNR_1[i] <- SNR(r)
  
}
#### SNR_2
alpha_2 <- 0.008+max_alpha[6]
SNR_2 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_2*r-L0+G-TN
}
r_2 <- seq(100,257,length.out = 106)
for (i in 1:length(r_2)) {
  r <- r_2[i]
  SNR_2[i] <- SNR(r)
  
}
#### SNR_2.75
alpha_2.75 <- 0.008+max_alpha[7]
SNR_2.75 <- vector("numeric")
SNR <- function(r){
  SNR <- SL-20*log10(r)-alpha_2.75*r-L0+G-TN
}
r_2.75 <- seq(100,200,length.out = 106)
for (i in 1:length(r_2.75)) {
  r <- r_2.75[i]
  SNR_2.75[i] <- SNR(r)
  
}

#### SNR_ ggplot for all distance for the scenod's peak
tiff("SNR_300_peak2.tiff", units="in", width=10, height=5, res=300)
plot(r_0,SNR_0,type="l",axes=F,ylim = c(0,50),xlim = c(100,2400),col="black",
     ann=F,xlab = NULL, ylab = NULL,lwd=2)
lines(r_0.01,SNR_0.01,type = "l",col="red",lwd=2)
lines(r_0.1,SNR_0.1,type = "l",col="blue",lwd=2)
lines(r_0.2,SNR_0.2,type = "l",col="green",lwd=2)
lines(r_0.5,SNR_0.5,type = "l",col="chocolate",lwd=2)
lines(r_1,SNR_1,type = "l",col="deeppink4",lwd=2)
lines(r_2,SNR_2,type = "l",col="#999900",lwd=2)
lines(r_2.75,SNR_2.75,type = "l",col="darkviolet",lwd=2)
points(100,SNR_12.67_peak1,col=434,lwd=2)
text(1900,18,"<SSC>=0",srt=-15,cex = 0.8,col = "black",lwd=2)
text(1450,18,"<SSC>=0.01(kg m-1)",srt=-14,cex = 0.8,col = "red",lwd=2)
text(1200,20,"<SSC>=0.1",srt=-20,cex = 0.8,col = "blue",lwd=2)
text(1000,18,"<SSC>=0.2",srt=-28,cex = 0.8,col = "green",lwd=2)
text(600,18,"<SSC>=0.5",srt=-50,cex = 0.8,col = "chocolate",lwd=2)
text(400,18,"<SSC>=1",srt=-60,cex = 0.8,col = "deeppink4",lwd=2)
text(250,18,"<SSC>=2",srt=-75,cex = 0.8,col = "#999900",lwd=2)
text(120,18,"<SSC>=2.75",srt=-75,cex = 0.8,col = "darkviolet",lwd=2)
text(85,18,"<SSC>=6.6",srt=-75,cex = 0.8,col = 434,lwd=2)
axis(side = 4, ylim=c(0,50),col="black",lwd=2,line = -1.5)
mtext(side = 4,text="SNR (dB)",line=0.5)
par(new=T)
plot(x_all,U_r, axes=F, ylim=c(0,max(U_r)), xlab="", ylab="",
     type="l",col="blue", main="",xlim=c(100,2400),lwd=2,lty = "dashed")
axis(side=2,ylim=c(0,max(U_r)),lwd=2,col = "blue",col.axis="blue",line = 1)
mtext(2,text = "Ur (m/s)",line=2.9,col = "blue")
axis(1,pretty(range(100,2400),10))
mtext("Range (m)",side=1,col="black",line=2)
dev.off()
### The maximum applicable measuremant distance (MAMD)...
###  Figure 5 of manuscript
SSC_1 <- c(12.67,5.28,2.75,2,1,0.5,0.2,0.1,0.01,0)
R_1 <- c(100,200,327,412,660.7,988,1473,1501,2260.5,2330)
SSC_2 <- c(6.67,2.75,2,1,0.5,0.2,0.1,0.01,0)
R_2 <- c(100,200,257,430,683,1140,1500,2200,2330)
tiff("Range.tiff", units="in", width=10, height=5, res=300)
plot(R_1,SSC_1,type = "l",lwd=2,ann=F,xlab = NULL, ylab = NULL,axes=F)
lines(R_2,SSC_2,lwd=2,col=2)
axis(1,pretty(range(100,2400),10))
mtext("Range (m)",side=1,col="black",line=2)
axis(side=2,lwd=2,col = "blue",line = 1,pretty(range(0,13),13),col.axis="blue")
mtext(2,text = "<SSC> (Kg m-1)",line=2.9,col = "blue")

dev.off()
