
packages <- c("refund", "dplyr", "refund.shiny")
## Now load or install&load all
package_check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
##### FPCA Example on real data #####

data(DTI)
MS <- subset(DTI, case ==1)  # subset data with multiple sclerosis (MS) case

index.na <- which(is.na(MS$cca))
Y <- MS$cca; Y[index.na] <- fpca.sc(Y)$Yhat[index.na]; sum(is.na(Y))
id <- MS$ID
visit.index <- MS$visit
visit.time <- MS$visit.time/max(MS$visit.time)

lfpca.dti1 <- fpca.lfda(Y = Y, subject.index = id,
                        visit.index = visit.index, obsT = visit.time,
                        LongiModel.method = 'lme',
                        mFPCA.pve = 0.95)
plot_shiny(lfpca.dti1)

lfpca.dti2 <- fpca.lfda(Y = Y, subject.index = id,
                        visit.index = visit.index, obsT = visit.time,
                        LongiModel.method = 'fpca.sc',
                        mFPCA.pve = 0.80, sFPCA.pve = 0.80)
plot_shiny(lfpca.dti2)

