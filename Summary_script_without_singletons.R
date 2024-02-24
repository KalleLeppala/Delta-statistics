#################Symetric####################
df1 <- read.csv(filename, header = TRUE, sep = ",", dec = ".")
df1 <- df1[, c(2, 4, 6, 8)]
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

for (i in seq(1, NROW(df1))) {
  for (j in seq(1, NCOL(df1))) {
    df1[i, j] <- signify(df1[i, j])
  }
}
#classifying each row <-> getting the number of windows
classify <- function(x) {
  y <- "unexplainable"
  if (prod(x == c( 0, 0, 0, 0))) {y <- "nothing"}
  if (prod(x == c( 1, 1, 1, 0))) {y <- "1 -> 3"}
  if (prod(x == c( 1, 0, 1, 1))) {y <- "3 -> 1"}
  if (prod(x == c(-1,-1, 0, 1))) {y <- "1 -> 4"}
  if (prod(x == c(-1, 0, 1, 1))) {y <- "4 -> 1"}
  if (prod(x == c( 1, 1,-1, 0))) {y <- "2 -> 3"}
  if (prod(x == c( 0, 1,-1,-1))) {y <- "3 -> 2"}
  if (prod(x == c(-1,-1, 0,-1))) {y <- "2 -> 4"}
  if (prod(x == c( 0,-1,-1,-1))) {y <- "4 -> 2"}
  if (prod(x == c( 0, 0,-1,-1))) {y <- "34 <-> 2*"} 
  if (prod(x == c( 0, 0, 1, 1))) {y <- "34 <-> 1*"}
  if (prod(x == c(-1,-1, 0, 0))) {y <- "12 <-> 4*"}
  if (prod(x == c( 1, 1, 0, 0))) {y <- "12 <->3*"}
  return(y)
}


# Create a vector -> store class labels for each row
class1 <- character(NROW(df1))

for (i in seq_along(class1 )) {
  class1 [i] <- classify(df1[i, ])
}

# Specify the list of scenarios -> include
specified_scenarios <- c("unexplainable", "nothing", "1 -> 3","3 -> 1","1 -> 4","4 -> 1","2 -> 3","2 -> 4","3 -> 2","4 -> 2","34 <-> 2*","34 <-> 1*","12 <-> 4*","12 <-> 3*")

# Create a table of class counts for specified scenarios
counts_class1 <- table(class1)

# Create a data frame with scenario labels <-> counts
table_class1 <- data.frame(Scenario = specified_scenarios, Count = 0)

# Update counts for existing scenarios using a loop
for (i in seq_along(specified_scenarios)) {
  scenario <- specified_scenarios[i]
  if (scenario %in% names(counts_class1)) {
    table_class1$Count[table_class1$Scenario == scenario] <- counts_class1[scenario]
  }
}

# Order the table by Count
filename2 <- paste("DS10mb_Summary_",combo, ".csv", sep = "")
write.csv(table_class1,filename2)

#################Asymetric####################
df1 <- read.csv(filename1, header = TRUE, sep = ",", dec = ".")
df1 <- df1[, c(2, 4, 6)]

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

for (i in seq(1, NROW(df1))) {
  for (j in seq(1, NCOL(df1))) {
    df1[i, j] <- signify(df1[i, j])
  }
}
#classifying each row and getting the number of windows
classify <- function(x) {
  y <- "unexplainable"
  if (prod(x == c( 0, 0, 0))) {y <- "nothing"}
  if (prod(x == c( 1, 1, 0))) {y <- "123 -> 2*"}
  if (prod(x == c(-1,-1, 0))) {y <- "123 -> 1*"}
  if (prod(x == c(-1, 0, 1))) {y <- "1 -> 4"}
  if (prod(x == c(-1,-1, 1))) {y <- "4 -> 1"}
  if (prod(x == c( 1, 0,-1))) {y <- "2 -> 4"}
  if (prod(x == c( 1, 1,-1))) {y <- "4 -> 2"}
  if (prod(x == c( 0,-1,-1))) {y <- "1 -> 5"}
  if (prod(x == c(-1,-1,-1))) {y <- "5 -> 1*"}
  if (prod(x == c( 0, 1, 1))) {y <- "2 -> 5"}
  if (prod(x == c( 1, 1, 1))) {y <- "5 -> 2*"}
  return(y)
}

# Create a vector -> store class labels for each row
class1 <- character(NROW(df1))

for (i in seq_along(class1)) {
  class1[i] <- classify(df1[i, ])
}

# Specify the list of scenarios -> include
specified_scenarios <- c("unexplainable", "nothing", "123 -> 2*","123 -> 1*","1 -> 4","4 -> 1","2 -> 4","4 -> 2","1 -> 5","5 -> 1*","2 -> 5","5 -> 2*")

# Create a table of class counts for specified scenarios
counts_class1 <- table(class1)

# Create a data frame with scenario labels and counts
table_class1 <- data.frame(Scenario = specified_scenarios, Count = 0)

# Update counts for existing scenarios using a loop
for (i in seq_along(specified_scenarios)) {
  scenario <- specified_scenarios[i]
  if (scenario %in% names(counts_class1)) {
    table_class1$Count[table_class1$Scenario == scenario] <- counts_class1[scenario]
  }
}

filename2 <- paste("DA5mb_Summary",combo, ".csv", sep = "")
write.csv(table_class1, filename2)

#################Quasiymetric####################
Q5MB_df1 <- read.csv(Q5MB_filename1, header = TRUE, sep = ",", dec = ".")
Q5MB_df1 <- Q5MB_df1[, c(2, 4, 6, 8)]

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

for (i in seq(1, NROW(Q5MB_df1))) {
  for (j in seq(1, NCOL(Q5MB_df1))) {
    Q5MB_df1[i, j] <- signify(Q5MB_df1[i, j])
  }
}
#classifying each row and getting the number of windows
classify <- function(x) {
  y <- "unexplainable"
  if (prod(x == c( 0, 0, 0, 0))) {y <- "nothing"}
  if (prod(x == c( 1, 1, 1, 0))) {y <- "1 -> 4"}
  if (prod(x == c( 1, 0, 1, 1))) {y <- "4 -> 1"}
  if (prod(x == c(-1,-1, 0, 1))) {y <- "1 -> 5"}
  if (prod(x == c(-1, 0, 1, 1))) {y <- "5 -> 1"}
  if (prod(x == c( 1, 1,-1, 0))) {y <- "2 -> 4"}
  if (prod(x == c( 0, 1,-1,-1))) {y <- "4 -> 2"}
  if (prod(x == c(-1,-1, 0,-1))) {y <- "2 -> 5"}
  if (prod(x == c( 0,-1,-1,-1))) {y <- "5 -> 2"}
  if (prod(x == c( 0, 0,-1,-1))) {y <- "45 <-> 2*"} 
  if (prod(x == c( 0, 0, 1, 1))) {y <- "45 <-> 1*"}
  if (prod(x == c(-1,-1, 0, 0))) {y <- "12 <-> 5*"}
  if (prod(x == c( 1, 1, 0, 0))) {y <- "12 <-> 4*"}
  return(y)
}

# Create a vector to store class labels for each row
Q5MB_class1 <- character(NROW(Q5MB_df1))

for (i in seq_along(Q5MB_class1)) {
  Q5MB_class1[i] <- classify(Q5MB_df1[i, ])
}

specified_scenarios <- c("unexplainable","nothing","1 -> 4","4 -> 1","1 -> 5","5 -> 1", "2 -> 4","4 -> 2","2 -> 5","5 -> 2","45 <-> 2*","45 <-> 1*","12 <-> 5*","12 <-> 4*","123 <-> 4*")

Q5MB_counts_class1 <- table(Q5MB_class1)

Q5MB_table_class1 <- data.frame(Scenario = specified_scenarios, Count = 0)

# Update counts for existing scenarios using a loop
for (i in seq_along(specified_scenarios)) {
  scenario <- specified_scenarios[i]
  if (scenario %in% names(Q5MB_counts_class1)) {
    Q5MB_table_class1$Count[Q5MB_table_class1$Scenario == scenario] <- Q5MB_counts_class1[scenario]
  }
}

# Order the table by Count
Q5MB_table_class1 <- Q5MB_table_class1[order(Q5MB_table_class1$Count, decreasing = TRUE), ]

Q5MB_filename2 <- paste("Summary_DQ500kb_new",combo, ".csv", sep = "")
write.csv(Q5MB_table_class1, Q5MB_filename2)
