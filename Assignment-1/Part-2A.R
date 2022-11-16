if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("yukai-yang/FE")
library(FE)
ls( grep("FE", search()) )

################################################################
##### Part 2: Asset Return Predictability and Market Efficiency
################################################################
## A. Testing for asset return predictability

acf(DJ_d$r_Dow_Jones[!is.na(DJ_d$r_Dow_Jones)],
    main="ACF of daily Dow Jones returns");
pacf(DJ_d$r_Dow_Jones[!is.na(DJ_d$r_Dow_Jones)],
     main="PACF of daily Dow Jones returns")

