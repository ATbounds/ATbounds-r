## code to prepare electronic fetal monitoring (EFM) and cesarean section (CS) Study dataset

rm(list = ls())

dat <- read.csv("./data-raw/efm.csv")

df <- {}

if (colnames(dat)[1] != 'count'){
  stop("The first column of the dataset must be 'count'.")   
}

# sample size
n <- sum(dat[,1])
if (n != 14484){
  stop("The sample size must be 14,484.")   
}

for(i in 1:nrow(dat)){
  count_i <- as.integer(dat[i,1])  # count of each cell
  vector_i <- as.matrix(dat[i,-1])  # type of each cell
  # repeat each cell type by its count
  df_i <- t(matrix(rep(vector_i,count_i), ncol = count_i, nrow = length(vector_i)))
  # combine observations
  df <- rbind(df, df_i)
}

df = df[sample.int(n, size = n),] # shuffle the rows 

colnames(df) <- colnames(dat)[-1]

EFM <- as.data.frame(df)

usethis::use_data(EFM, overwrite = TRUE)