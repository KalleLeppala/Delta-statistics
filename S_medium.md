Smaller simulation in the symmetric tree $S$
================
Kalle Leppälä

This is a smaller version of [S_big.md](S_big.md) with the number of
independent allelic patterns reduced to 100,000.

## Simulating the 32 admixture scenarios in the symmetric tree $S$

``` r
tables <- list()
lambda <- 1
mu <- 0.0001
N <- 100000

# 1 -> 3
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1), LL = c(0.5, 0.1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 3 -> 1
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1), RL = c(0.5, 0.1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 1 -> 4
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1), LL = c(0.5, 0.1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 4 -> 1
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1), RR = c(0.5, 0.1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 2 -> 3
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1), LR = c(0.5, 0.1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 3 -> 2
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1), RL = c(0.5, 0.1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 2 -> 4
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1), LR = c(0.5, 0.1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 4 -> 2
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RL = c(2, 1), RR = c(0.5, 0.1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(OG = c(3, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(RL = c(1, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 12 -> 3
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(0.75, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2.25, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1.5, 1), LR = c(1.125, 0.1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 3 -> 12
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(0.75, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2.25, 1), RL = c(1.125, 0.1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1.5, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 12 -> 4
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(0.75, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2.25, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1.5, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1), LR = c(1.125, 0.1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 4 -> 12
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(0.75, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RL = c(2.25, 1), RR = c(1.125, 0.1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(OG = c(3, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(RL = c(1.5, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 34 -> 1
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1.5, 1), RR = c(1.125, 0.1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2.25, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(0.75, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 1 -> 34
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1.5, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2.25, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(0.75, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1), LL = c(1.125, 0.1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 34 -> 2
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(RR = c(2.25, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(LL = c(1.5, 1), RR = c(1.125, 0.1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(0.75, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 2 -> 34
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1.5, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2.25, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(0.75, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1), LR = c(1.125, 0.1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 1 -> 5
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list(LL = c(0.5, 0.1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 5 -> 1
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1), OG = c(0.5, 0.1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 2 -> 5
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list(LR = c(0.5, 0.1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 5 -> 2
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1), OG = c(0.5, 0.1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 3 -> 5
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list(RL = c(0.5, 0.1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 5 -> 3
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1), OG = c(0.5, 0.1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 4 -> 5
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list(RR = c(0.5, 0.1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 5 -> 4
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RL = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(OG = c(3, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(RL = c(1, 1), OG = c(0.5, 0.1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 12345 -> 1
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1), GH = c(0.5, 0.1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
GH <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(), switches = list(OG = c(3.5, 1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG, GH = GH)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 12345 -> 2
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1), GH = c(0.5, 0.1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
GH <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(), switches = list(OG = c(3.5, 1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG, GH = GH)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 12345 -> 3
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1), GH = c(0.5, 0.1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
GH <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(), switches = list(OG = c(3.5, 1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG, GH = GH)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 12345 -> 4
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1), GH = c(0.5, 0.1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
GH <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(), switches = list(OG = c(3.5, 1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG, GH = GH)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 1234 -> 1
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1), GH = c(0.5, 0.1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
GH <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(), switches = list(RR = c(2.5, 1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG, GH = GH)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 1234 -> 2
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1), GH = c(0.5, 0.1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
GH <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(), switches = list(RR = c(2.5, 1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG, GH = GH)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 1234 -> 3
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1), GH = c(0.5, 0.1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
GH <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(), switches = list(RR = c(2.5, 1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG, GH = GH)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table
# 1234 -> 4
LL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P1 = 0), switches = list(LR = c(1, 1)))
LR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P2 = 0), switches = list(RR = c(2, 1)))
RL <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P3 = 0), switches = list(RR = c(1, 1)))
RR <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P4 = 0), switches = list(OG = c(3, 1), GH = c(0.5, 0.1)))
OG <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(P5 = 0), switches = list())
GH <- list(lambda = lambda, mu = mu, inhabitants = character(0), samples = list(), switches = list(RR = c(2.5, 1)))
graph <- list(LL = LL, LR = LR, RL = RL, RR = RR, OG = OG, GH = GH)
abbababa <- character(0)
for (k in seq(1, N)) {
  abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)
    while (is.na(abbababa[k]) == TRUE) {abbababa[k] <- pattern(solve_tree(graph), Tmax = 0.0055)}
}
table <- table(abbababa)
tables[[length(tables) + 1]] <- table

save(tables, file = "S_medium.RData")
```

## Computing the $\Delta$-statistics

``` r
load(file = "S_medium.RData")
event <- character(0)
statistic <- character(0)
value <- numeric(0)
zvalue <- numeric(0)
sign <- numeric(0)
ci <- numeric(0)

# 1 -> 3
table <- tables[[1]]
event <- c(event, rep("1 \u2794 3", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 3 -> 1
table <- tables[[2]]
event <- c(event, rep("3 \u2794 1", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 1 -> 4
table <- tables[[3]]
event <- c(event, rep("1 \u2794 4", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 4 -> 1
table <- tables[[4]]
event <- c(event, rep("4 \u2794 1", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 2 -> 3
table <- tables[[5]]
event <- c(event, rep("2 \u2794 3", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 3 -> 2
table <- tables[[6]]
event <- c(event, rep("3 \u2794 2", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 2 -> 4
table <- tables[[7]]
event <- c(event, rep("2 \u2794 4", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 4 -> 2
table <- tables[[8]]
event <- c(event, rep("4 \u2794 2", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 12 -> 3
table <- tables[[9]]
event <- c(event, rep("12 \u2794 3", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 3 -> 12
table <- tables[[10]]
event <- c(event, rep("3 \u2794 12", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 12 -> 4
table <- tables[[11]]
event <- c(event, rep("12 \u2794 4", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 4 -> 12
table <- tables[[12]]
event <- c(event, rep("4 \u2794 12", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 34 -> 1
table <- tables[[13]]
event <- c(event, rep("34 \u2794 1", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 1 -> 34
table <- tables[[14]]
event <- c(event, rep("1 \u2794 34", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 34 -> 2
table <- tables[[15]]
event <- c(event, rep("34 \u2794 2", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 2 -> 34
table <- tables[[16]]
event <- c(event, rep("2 \u2794 34", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 1 -> 5
table <- tables[[17]]
event <- c(event, rep("1 \u2794 5", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 5 -> 1
table <- tables[[18]]
event <- c(event, rep("5 \u2794 1", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 2 -> 5
table <- tables[[19]]
event <- c(event, rep("2 \u2794 5", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 5 -> 2
table <- tables[[20]]
event <- c(event, rep("5 \u2794 2", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 3 -> 5
table <- tables[[21]]
event <- c(event, rep("3 \u2794 5", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 5 -> 3
table <- tables[[22]]
event <- c(event, rep("5 \u2794 3", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 4 -> 5
table <- tables[[23]]
event <- c(event, rep("4 \u2794 5", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 5 -> 4
table <- tables[[24]]
event <- c(event, rep("5 \u2794 4", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 12345 -> 1
table <- tables[[25]]
event <- c(event, rep("12345 \u2794 1", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 12345 -> 2
table <- tables[[26]]
event <- c(event, rep("12345 \u2794 2", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 12345 -> 3
table <- tables[[27]]
event <- c(event, rep("12345 \u2794 3", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 12345 -> 4
table <- tables[[28]]
event <- c(event, rep("12345 \u2794 4", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 1234 -> 1
table <- tables[[29]]
event <- c(event, rep("1234 \u2794 1", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 1234 -> 2
table <- tables[[30]]
event <- c(event, rep("1234 \u2794 2", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 1234 -> 3
table <- tables[[31]]
event <- c(event, rep("1234 \u2794 3", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])
# 1234 -> 4
table <- tables[[32]]
event <- c(event, rep("1234 \u2794 4", 8))
statistic <- c(statistic, "DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268")
value <- c(value, DS16(table)[[1]], DS26(table)[[1]], DS35(table)[[1]], DS45(table)[[1]],
                    DS57(table)[[1]], DS68(table)[[1]], DS3457(table)[[1]], DS1268(table)[[1]])
zvalue <- c(zvalue, DS16(table)[[2]], DS26(table)[[2]], DS35(table)[[2]], DS45(table)[[2]],
                      DS57(table)[[2]], DS68(table)[[2]], DS3457(table)[[2]], DS1268(table)[[2]])
sign <- c(sign, DS16(table)[[3]], DS26(table)[[3]], DS35(table)[[3]], DS45(table)[[3]],
                  DS57(table)[[3]], DS68(table)[[3]], DS3457(table)[[3]], DS1268(table)[[3]])
ci <- c(ci, DS16(table)[[4]], DS26(table)[[4]], DS35(table)[[4]], DS45(table)[[4]],
              DS57(table)[[4]], DS68(table)[[4]], DS3457(table)[[4]], DS1268(table)[[4]])

symmetric <- data.frame(event = event, statistic = statistic, value = value, zvalue = zvalue, sign = sign, ci = ci)
symmetric$statistic <- factor(symmetric$statistic, levels <- c("DS16", "DS26", "DS35", "DS45", "DS57", "DS68", "DS3457", "DS1268"))
symmetric$sign <- sub("-", "\u2212", symmetric$sign) # A hyphen is not a minus sign.
save(symmetric, file = "S_medium_stats.RData")
```

## Plotting the results

``` r
load(file = "S_medium_stats.RData")
purples <- c(alpha("#810F7C", 1),
             alpha("#810F7C", 0.9),
             alpha("#810F7C", 0.8),
             alpha("#810F7C", 0.7),
             alpha("#810F7C", 0.6),
             alpha("#810F7C", 0.5),
             alpha("#810F7C", 0.4),
             alpha("#810F7C", 0.3))
labels <- c(expression(paste(""[italic(S)],Delta,"*"[1-6])),
            expression(paste(""[italic(S)],Delta,"*"[2-6])),
            expression(paste(""[italic(S)],Delta,"*"[3-5])),
            expression(paste(""[italic(S)],Delta,"*"[4-5])),
            expression(paste(""[italic(S)],Delta,"*"[5+7])),
            expression(paste(""[italic(S)],Delta,"*"[6+8])),
            expression(paste(""[italic(S)],Delta,"*"[3+4-5+7])),
            expression(paste(""[italic(S)],Delta,"*"[1+2-6+8])))

events <- c("1 \u2794 3", "3 \u2794 1", "1 \u2794 4", "4 \u2794 1",
            "2 \u2794 3", "3 \u2794 2", "2 \u2794 4", "4 \u2794 2")
temp <- symmetric[symmetric$event %in% events, ]
temp$event <- factor(temp$event, levels <- events)
plot1 <- ggplot(temp, aes(fill = statistic, y = value, x = event)) +
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  theme(legend.text.align = 0, legend.key.size = unit(0.22, 'cm')) +
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = c(-0.5, 0, 0.5), labels = c("\u22120.5", "0.0", "0.5")) +
  scale_fill_manual(values = purples, name = "", labels = labels) +
  geom_hline(yintercept = 0) +
    geom_bar(width = 0.8, position = position_dodge(0.8), stat = "identity") +
    geom_text(position = position_dodge(0.8), aes(y = 0.2, label = sign, hjust = 0.5), size = 2.8)
```

    ## Warning: The `legend.text.align` argument of `theme()` is deprecated as of ggplot2
    ## 3.5.0.
    ## ℹ Please use theme(legend.text = element_text(hjust)) instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
leg <- as_ggplot(get_legend(plot1))
plot1 <- plot1 + theme(legend.position = "none")

events <- c("1 \u2794 5", "5 \u2794 1", "2 \u2794 5", "5 \u2794 2",
            "3 \u2794 5", "5 \u2794 3", "4 \u2794 5", "5 \u2794 4")
temp <- symmetric[symmetric$event %in% events, ]
temp$event <- factor(temp$event, levels <- events)
plot2 <- ggplot(temp, aes(fill = statistic, y = value, x = event)) +
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = c(-0.5, 0, 0.5), labels = c("\u22120.5", "0.0", "0.5")) +
  scale_fill_manual(values = purples, name = "", labels = labels) +
  geom_hline(yintercept = 0) +
    geom_bar(width = 0.8, position = position_dodge(0.8), stat = "identity") +
    geom_text(position = position_dodge(0.8), aes(y = 0.2, label = sign, hjust = 0.5), size = 2.8)

events <- c("1234 \u2794 1", "1234 \u2794 2", "1234 \u2794 3", "1234 \u2794 4",
            "12345 \u2794 1", "12345 \u2794 2", "12345 \u2794 3", "12345 \u2794 4")
temp <- symmetric[symmetric$event %in% events, ]
temp$event <- factor(temp$event, levels <- events)
plot3 <- ggplot(temp, aes(fill = statistic, y = value, x = event)) +
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = c(-0.5, 0, 0.5), labels = c("\u22120.5", "0.0", "0.5")) +
  scale_fill_manual(values = purples, name = "", labels = labels) +
  geom_hline(yintercept = 0) +
    geom_bar(width = 0.8, position = position_dodge(0.8), stat = "identity") +
    geom_text(position = position_dodge(0.8), aes(y = 0.2, label = sign, hjust = 0.5), size = 2.8)

events <- c("12 \u2794 3", "3 \u2794 12", "12 \u2794 4", "4 \u2794 12",
            "34 \u2794 1", "1 \u2794 34", "34 \u2794 2", "2 \u2794 34")
temp <- symmetric[symmetric$event %in% events, ]
temp$event <- factor(temp$event, levels <- events)
plot4 <- ggplot(temp, aes(fill = statistic, y = value, x = event)) +
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = c(-0.5, 0, 0.5), labels = c("\u22120.5", "0.0", "0.5")) +
  scale_fill_manual(values = purples, name = "", labels = labels) +
  geom_hline(yintercept = 0) +
    geom_bar(width = 0.8, position = position_dodge(0.8), stat = "identity") +
    geom_text(position = position_dodge(0.8), aes(y = 0.2, label = sign, hjust = 0.5), size = 2.8)

figure <- ggarrange(plot1, plot2, plot3, plot4, nrow = 4, ncol = 1)
figure <- ggarrange(figure, leg, nrow = 1, ncol = 2, widths = c(7, 1))
ggsave("S_medium.png", plot = figure, width = 160, height = 80, units = "mm")
print(figure)
```

![](S_medium_files/figure-gfm/plots-1.png)<!-- -->
