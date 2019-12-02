#### load packages
library(flowCore)
library(flowViz)
library(gdata)
library(Hmisc)

xvalues <- c(1, 1.75, 3.5, 7, 15, 30, 62, 125, 250, 500, 1000, 2000)
conc = as.character (c(1, 1.75, 3.5, 7, 15, 30, 62, 125, 250, 500, 1000, 2000))

FRY418t24h <- read.flowSet(files=NULL, "FC_LEXA-ER-B42/FC_LEXA-ER-B42-raw")

# define the covariance matrix with the calculated values
cov <- matrix(c(5913616.473, 1173108.771, 1173108.771, 2158883.527), ncol=2, dimnames=list(c("FSC-W", "SSC-W"), c("FSC-W", "SSC-W")))
# define center of elipsoid
center <- c("FSC-W"=67200, "SSC-W"=66200) 
# gate
elliGate <- ellipsoidGate(filterId= "ElliGate", .gate=cov, mean=center)

##### data extraction

# FRY418
gFRY418t24h <- filter(FRY418t24h, elliGate)
# store gated samples as new frame
gatedFRY418t24h <- Subset(FRY418t24h, gFRY418t24h)

#### fluorescence analysis
FRY418t24hmKate2 <- as.data.frame(gatedFRY418t24h[["1.fcs"]]@exprs[,"561 [B]-A"])
FRY418t24hcitrine <- as.data.frame(gatedFRY418t24h[["1.fcs"]]@exprs[,"488 [F]-A"])

for (i in 2:12){
  # mKate2
  xx <- as.data.frame(gatedFRY418t24h[[paste(as.character(i),".fcs", sep="")]]@exprs[,"561 [B]-A"])
  FRY418t24hmKate2 <- cbindX(FRY418t24hmKate2,xx)
  # citrine
  xx <- as.data.frame(gatedFRY418t24h[[paste(as.character(i),".fcs", sep="")]]@exprs[,"488 [F]-A"])
  FRY418t24hcitrine <- cbindX(FRY418t24hcitrine,xx)
}


colnames(FRY418t24hmKate2)<- paste("FRY418t24hmKate2.",conc,sep="")
attach(FRY418t24hmKate2)
colnames(FRY418t24hcitrine)<- paste("FRY418t24hcitrine.",conc,sep="")
attach(FRY418t24hcitrine)

FRY418t24hmedianmKate2<- apply(FRY418t24hmKate2, 2, median, na.rm=TRUE)
FRY418t24hq25mKate2 <- apply(FRY418t24hmKate2, 2, quantile, prob=0.25, na.rm=TRUE)
FRY418t24hq75mKate2 <- apply(FRY418t24hmKate2, 2, quantile, prob=0.75, na.rm=TRUE)

FRY418t24hmedianCitrine<- apply(FRY418t24hcitrine, 2, median, na.rm=TRUE)
FRY418t24hq25Citrine <- apply(FRY418t24hcitrine, 2, quantile, prob=0.25, na.rm=TRUE)
FRY418t24hq75Citrine <- apply(FRY418t24hcitrine, 2, quantile, prob=0.75, na.rm=TRUE)






########################################## plots
# Citrine
pdf(file="FC_LEXA-ER-B42-analyzed/FRY418t24hCitrine.pdf")
par(fg="black")
plot(x=xvalues, 
     y=FRY418t24hmedianCitrine, 
     ylim=c(2,60000),
     log="x",
     main="FRY418 24 h induction",
     xlab="[beta-estradiol] (nM)", 
     ylab="Citrine fluorescences (AU)", 
     col="blue", pch=18)
par(fg="blue")
errbar(x=xvalues, y=FRY418t24hmedianCitrine, yplus=FRY418t24hq25Citrine, yminus=FRY418t24hq75Citrine, add=TRUE, col="blue", pch=18)

dev.off()


# mKate2
pdf(file="FC_LEXA-ER-B42-analyzed/FRY418t24hmKate2.pdf")
par(fg="black")
plot(x=xvalues, 
     y=FRY418t24hmedianmKate2, 
     ylim=c(2,7000),
     log="x",
     main="FRY418 24 h induction",
     xlab="[beta-estradiol] (nM)", 
     ylab="mKate2 fluorescences (AU)", 
     col="blue", pch=18)
par(fg="blue")
errbar(x=xvalues, y=FRY418t24hmedianmKate2, yplus=FRY418t24hq25mKate2, yminus=FRY418t24hq75mKate2, add=TRUE, col="blue", pch=18)

dev.off()
