library(FE)


# 1. ####

# Here we transform portcap_m to label dummies
# for quantiles of marketcap-weighted returns

# make data general for function wrapper
data <- portcap_m

cols <- colnames(data)

# keep only segmented returns columns
cols <- cols[2: (length(cols)-2) ]

n <- nrow(data)
m <- length(cols)

res <- matrix(0,
              nrow=n*m,
              ncol=m+2,
)

# set column names
colnames(res) <- c('Date', 'r', cols)

# initiate returns with NAs
res[,'r'] <- NA

# iterate for each column
for (i in 1:m) {

  col_data <- t(portcap_m[, col])

  start <- 1 + ((i-1)*n)
  end <- i * n
  col <- cols[i]

  # set return
  res[start:end, 'r'] <- col_data

  # set date
  res[start:end, 'Date'] <- t(portcap_m[, 'Date'])

  # set dummy to True
  res[start:end, col] <- 1

}

# save as dataframe
df <- data.frame(res)

quantiles <- list(
        1:3, # 30% 40% 30% quantiles
        4:8, # 20% quantiles
        9:18, # 10% deciles
        1:length(cols) # all
        )

colnames(data)

for (q in quantiles){

  # make lm formula
  dep_vars <- paste(cols[q], collapse=' + ')
  formula <- paste('r ~ ', dep_vars)

  # estimate regression
  model <- lm(formula=formula, data=df)
  print(summary(model))

}

view(df)
