#################Symetric####################
S500kb_df1 <- read.csv("DS500kb_pABC_BC,pEBB,pAK,pPB,BLK.csv", header = TRUE, sep = ",", dec = ".")
S500kb_df1 <- S500kb_df1[, c(2, 4, 6, 8, 10, 12, 14, 16)]
#converting values
signify <- function(x) {
  y <- 0
  if (x == "'NA'") {y <- -9}
  else {
    if (as.numeric(x) > 2.575829) {y <- 1}
    if (as.numeric(x) < -2.575829) {y <- -1}
  }
  return(y)
}

for (i in seq(1, NROW(S500kb_df1))) {
  for (j in seq(1, NCOL(S500kb_df1))) {
    S500kb_df1[i, j] <- signify(S500kb_df1[i, j])
  }
}
#classifying each row <-> getting the number of windows
classify <- function(x) {
  y <- "unexplainable"
  if (prod(x == c( 0, 0, 0, 0, 0, 0, 0, 0))) {y <- "nothing"}
  if (prod(x == c( 1, 1, 1, 0,-1,-1, 0, 0))) {y <- "1 -> 3"}
  if (prod(x == c( 1, 0, 1, 1,-1,-1, 0, 0))) {y <- "3 -> 1"}
  if (prod(x == c(-1,-1, 0, 1,-1, 1, 0, 0))) {y <- "1 -> 4"}
  if (prod(x == c(-1, 0, 1, 1,-1, 1, 0, 0))) {y <- "4 -> 1"}
  if (prod(x == c( 1, 1,-1, 0, 1,-1, 0, 0))) {y <- "2 -> 3"}
  if (prod(x == c( 0, 1,-1,-1, 1,-1, 0 ,0))) {y <- "3 -> 2"}
  if (prod(x == c(-1,-1, 0,-1, 1, 1, 0, 0))) {y <- "2 -> 4"}
  if (prod(x == c( 0,-1,-1,-1, 1, 1, 0, 0))) {y <- "4 -> 2"}
  if (prod(x == c( 0, 0,-1,-1, 1, 0, 1 ,0))) {y <- "12345-> 1"} 
  if (prod(x == c( 0, 0,-1,-1, 1, 0, 0, 0))) {y <- "34 <-> 2*"}
  if (prod(x == c( 0, 0,-1,-1, 1, 0,-1, 0))) {y <- "5 -> 1"}
  if (prod(x == c( 0, 0,-1,-1, 0, 0,-1, 0))) {y <- "1 -> 5"}
  if (prod(x == c( 0, 0, 1, 1,-1, 0,-1 ,0))) {y <- "12345 -> 2"} 
  if (prod(x == c( 0, 0, 1, 1,-1, 0, 0, 0))) {y <- "34 <-> 1*"}
  if (prod(x == c( 0, 0, 1, 1,-1, 0, 1, 0))) {y <- "5 -> 2"}
  if (prod(x == c( 0, 0, 1, 1, 0, 0, 1, 0))) {y <- "2 -> 5"}
  if (prod(x == c(-1,-1, 0, 0, 0, 1, 0, 1))) {y <- "12345 -> 3"} 
  if (prod(x == c(-1,-1, 0, 0, 0, 1, 0, 0))) {y <- "12 <-> 4*"}
  if (prod(x == c(-1,-1, 0, 0, 0, 1, 0,-1))) {y <- "5 -> 3"}
  if (prod(x == c(-1,-1, 0, 0, 0, 0, 0,-1))) {y <- "3 -> 5"}
  if (prod(x == c( 1, 1, 0, 0, 0,-1, 0,-1))) {y <- "12345 -> 4"} 
  if (prod(x == c( 1, 1, 0, 0, 0,-1, 0, 0))) {y <- "12 <-> 3*"}
  if (prod(x == c( 1, 1, 0, 0, 0,-1, 0, 1))) {y <- "5 -> 4"}
  if (prod(x == c( 1, 1, 0, 0, 0, 0, 0, 1))) {y <- "4 -> 5"}
  return(y)
}


# Create a vector -> store class labels for each row
S500kb_class1 <- character(NROW(S500kb_df1))

for (i in seq_along(S500kb_class1 )) {
  S500kb_class1 [i] <- classify(S500kb_df1[i, ])
}

# Specify the list of scenarios -> include
specified_scenarios <- c("unexplainable", "nothing", "1 -> 3","3 -> 1","1 -> 4","4 -> 1","2 -> 3","3 -> 2","2 -> 4","4 -> 2","12345 -> 1","34 <-> 2*","5 -> 1","1 -> 5","12345 -> 2","34 <-> 1*","5 -> 2","2 -> 5","12345 -> 3","12 <-> 4*","5 -> 3","3 -> 5","12345 -> 4","12 <-> 3*","5 -> 4","4 -> 5")

# Create a table of class counts for specified scenarios
S500kb_counts_class1 <- table(S500kb_class1)

# Create a data frame with scenario labels <-> counts
S500kb_table_class1 <- data.frame(Scenario = specified_scenarios, Count = 0)

# Update counts for existing scenarios using a loop
for (i in seq_along(specified_scenarios)) {
  scenario <- specified_scenarios[i]
  if (scenario %in% names(S500kb_counts_class1)) {
    S500kb_table_class1$Count[S500kb_table_class1$Scenario == scenario] <- S500kb_counts_class1[scenario]
  }
}

S500kb_filename2 <- paste("DS500KB_AK_", ".csv", sep = "")
write.csv(S500kb_table_class1, S500kb_filename2)

#################Asymetric####################
A500kb_df1 <- read.csv("DApABC_BC,pBB,pEBB,pPB,BLK.500kb.csv", header = TRUE, sep = ",", dec = ".")
A500kb_df1 <- A500kb_df1[, c(2, 4, 6, 8)]

#converting values
signify <- function(x) {
  y <- 0
  if (x == "NA'") {y <- -9}
  else {
    if (as.numeric(x) > 2.575829) {y <- 1}
    if (as.numeric(x) < -2.575829) {y <- -1}
  }
  return(y)
}

for (i in seq(1, NROW(A500kb_df1))) {
  for (j in seq(1, NCOL(A500kb_df1))) {
    A500kb_df1[i, j] <- signify(A500kb_df1[i, j])
  }
}
#classifying each row <-> getting the number of windows
classify <- function(x) {
  y <- "unexplainable"
  if (prod(x == c( 0, 0, 0, 0))) {y <- "nothing"}
  if (prod(x == c( 1, 1, 0, 0))) {y <- "1 <-> 3*"}
  if (prod(x == c(-1,-1, 0, 0))) {y <- "2 <-> 3*"}
  if (prod(x == c(-1, 0, 1, 0))) {y <- "1 -> 4"}
  if (prod(x == c(-1,-1, 1, 0))) {y <- "4 -> 1"}
  if (prod(x == c( 1, 0,-1, 0))) {y <- "2 -> 4"}
  if (prod(x == c( 1, 1,-1, 0))) {y <- "4 -> 2"}
  if (prod(x == c( 0,-1,-1,-1))) {y <- "1 -> 5"}
  if (prod(x == c(-1,-1,-1,-1))) {y <- "5 -> 1"}
  if (prod(x == c(-1,-1,-1, 0))) {y <- "1234 -> 1"} 
  if (prod(x == c(-1,-1,-1, 1))) {y <- "12345 -> 1"}
  if (prod(x == c( 0, 1, 1, 1))) {y <- "2 -> 5"}
  if (prod(x == c( 1, 1, 1, 1))) {y <- "5 -> 2"}
  if (prod(x == c( 1, 1, 1, 0))) {y <- "1234 -> 2"}
  if (prod(x == c( 1, 1, 1,-1))) {y <- "12345 -> 2"}
  return(y)
}

# Create a vector -> store class labels for each row
A500kb_class1 <- character(NROW(A500kb_df1))

for (i in seq_along(A500kb_class1)) {
  A500kb_class1[i] <- classify(A500kb_df1[i, ])
}

# Specify the list of scenarios -> include
specified_scenarios <- c("unexplainable", "nothing", "1 <-> 3*","2 <-> 3*","1 -> 4","4 -> 1","2 -> 4","4 -> 2","1 -> 5","5 -> 1","1234-> 1","12345 -> 1","2 -> 5", "5 -> 2", "1234 -> 2", "12345 -> 2")
# Create a table of class counts for specified scenarios
A500kb_counts_class1 <- table(A500kb_class1)

# Create a data frame with scenario labels <-> counts
A500kb_table_class1 <- data.frame(Scenario = specified_scenarios, Count = 0)

# Update counts for existing scenarios using a loop
for (i in seq_along(specified_scenarios)) {
  scenario <- specified_scenarios[i]
  if (scenario %in% names(A500kb_counts_class1)) {
    A500kb_table_class1$Count[A500kb_table_class1$Scenario == scenario] <- A500kb_counts_class1[scenario]
  }
}

A500kb_filename2 <- paste("DA500kb_plot_",combo, ".csv", sep = "")
write.csv(A500kb_table_class1, A500kb_filename2)


#################Quasisymetric####################
Q5mb <- read.csv("DQ5mb_pABC_BC,pBB,pEBB,pAK,pPB.csv", header = TRUE, sep = ",", dec = ".")
Q5mb_df1 <- data.frame(Q5mb)
Q5mb_df1 <- Q5mb_df1[, c(2, 4, 6, 8, 10, 12, 14, 16)]

#converting values
signify <- function(x) {
  y <- 0
  if (x == "NA'") {y <- -9}
  else {
    if (as.numeric(x) > 2.575829) {y <- 1}
    if (as.numeric(x) < -2.575829) {y <- -1}
  }
  return(y)
}

for (i in seq(1, NROW(Q5mb_df1))) {
  for (j in seq(1, NCOL(Q5mb_df1))) {
    Q5mb_df1[i, j] <- signify(Q5mb_df1[i, j])
  }
}
#classifying each row <-> getting the number of windows
classify <- function(x) {
  y <- "unexplainable"
  if (prod(x == c( 0, 0, 0, 0, 0, 0, 0, 0))) {y <- "nothing"}
  if (prod(x == c( 1, 1, 1, 0,-1,-1, 0,-1))) {y <- "1 -> 4"}
  if (prod(x == c( 1, 0, 1, 1,-1,-1, 1 ,0))) {y <- "4 -> 1"}
  if (prod(x == c(-1,-1, 0, 1,-1, 1, 0, 1))) {y <- "1 -> 5"}
  if (prod(x == c(-1, 0, 1, 1,-1, 1, 1, 0))) {y <- "5 -> 1"}
  if (prod(x == c( 1, 1,-1, 0, 1,-1, 0,-1))) {y <- "2 -> 4"}
  if (prod(x == c( 0, 1,-1,-1, 1,-1,-1, 0))) {y <- "4 -> 2"}
  if (prod(x == c(-1,-1, 0,-1, 1, 1, 0, 1))) {y <- "2 -> 5"}
  if (prod(x == c( 0,-1,-1,-1, 1, 1,-1, 0))) {y <- "5 -> 2"}
  if (prod(x == c (0, 0,-1,-1, 1, 0, 0, 0))) {y <- "2 -> 45"}
  if (prod(x == c( 0, 0,-1,-1, 1, 0,-1, 0))) {y <- "45 -> 2"}
  if (prod(x == c( 0, 0,-1,-1, 0, 0,-1, 0))) {y <- "1 <-> 3*"}
  if (prod(x == c( 0, 0,-1,-1,-1, 0,-1, 0))) {y <- "12345 -> 2"}
  if (prod(x == c( 0, 0, 1, 1,-1, 0, 0, 0))) {y <- "1 -> 45"}
  if (prod(x == c( 0, 0, 1, 1,-1, 0, 1, 0))) {y <- "45 -> 1"}
  if (prod(x == c( 0, 0, 1, 1, 0, 0, 1, 0))) {y <- "2 <-> 3*"}
  if (prod(x == c( 0, 0, 1, 1, 1, 0, 1, 0))) {y <- "12345 -> 1"}
  if (prod(x == c(-1,-1, 0, 0, 0, 1, 0, 0))) {y <- "5 -> 12"}
  if (prod(x == c(-1,-1, 0, 0, 0, 1, 0, 1))) {y <- "12 -> 5"}
  if (prod(x == c(-1,-1, 0, 0, 0, 0, 0,-1))) {y <- "4 -> 3"}
  if (prod(x == c(-1,-1, 0, 0, 0,-1, 0,-1))) {y <- "3 -> 4"}
  if (prod(x == c( 1, 1, 0, 0, 0,-1, 0, 0))) {y <- "4 -> 12"}
  if (prod(x == c( 1, 1, 0, 0, 0,-1, 0,-1))) {y <- "12 -> 4"}
  if (prod(x == c( 1, 1, 0, 0, 0, 0, 0, 1))) {y <- "5 -> 3"}
  if (prod(x == c( 1, 1, 0, 0, 0, 1, 0, 1))) {y <- "3 -> 5"}
  if (prod(x == c( 0, 0, 0, 0, 0, 1, 0, 1))) {y <- "123 <-> 5*"}
  if (prod(x == c( 0, 0, 0, 0, 0,-1, 0,-1))) {y <- "123 <-> 4*"}
  return(y)
}

# Create a vector -> store class labels for each row
Q5mb_class1 <- character(NROW(Q5mb_df1))

for (i in seq_along(Q5mb_class1)) {
  Q5mb_class1[i] <- classify(Q5mb_df1[i, ])
}

specified_scenarios <- c("unexplainable","nothing","1 -> 4","4 -> 1","1 -> 5","5 -> 1","2 -> 4","4 -> 2","2 -> 5","5 -> 2","2 -> 45","45 -> 2","1 <-> 3*","12345 -> 2","1 -> 45","45 -> 1","2 <-> 3*","12345 -> 1","5 -> 12","12 -> 5","4 -> 3","3 -> 4" ,"4 -> 12","12 -> 4","5 -> 3","3 -> 5","123 -> 5*","123 -> 4*")

Q5mb_counts_class1 <- table(Q5mb_class1)

Q5mb_table_class1 <- data.frame(Scenario = specified_scenarios, Count = 0)

# Update counts for existing scenarios using a loop
for (i in seq_along(specified_scenarios)) {
  scenario <- specified_scenarios[i]
  if (scenario %in% names(Q5mb_counts_class1)) {
    Q5mb_table_class1$Count[Q5mb_table_class1$Scenario == scenario] <- Q5mb_counts_class1[scenario]
  }
}

# Order the table by Count


library(dplyr)
Q5mb_filename2 <- paste("Summary_DQ5mb_new",combo, ".csv", sep = "")
write.csv(Q5mb_table_class1, Q5mb_filename2)

write.csv(df_class1, "pABC_A,pBB,pAPB,pPB,BLK_summary.csv")


