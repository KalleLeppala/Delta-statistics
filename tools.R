## ----setup, include = FALSE---------------------------------------------------

knitr::opts_chunk$set(echo = TRUE)
set.seed(0)



## ----echo = FALSE, out.width = "50%", fig.align = "center"--------------------

knitr::include_graphics("sampling.png")



## ----simulation functions-----------------------------------------------------

# Better documentation of these functions is found in the markdown document tools.md.

solve_branch <- function(branch, gene_tree, start, end) {
  # If there's at most one lineage in the branch during the time period, we only record the time passing.
  if (length(branch$inhabitants) <= 1) {
    for (inhabitant in branch$inhabitants) {
      if (inhabitant %in% names(gene_tree)) {gene_tree[[inhabitant]] <- gene_tree[[inhabitant]] + branch$mu*(end - start)}
      else {gene_tree[[inhabitant]] <- branch$mu*(end - start)}
    }
  }
  # If there's two or more lineages, coalescence events might occur.
  if (length(branch$inhabitants) > 1) {
    minimum_time <- rexp(1, rate = choose(length(branch$inhabitants), 2)*branch$lambda)
    if (minimum_time < end - start) {# If there's at least one coalescence event during the time frame, we record it, record the time passing, and call the function recursively.
      for (inhabitant in branch$inhabitants) {
        if (inhabitant %in% names(gene_tree)) {gene_tree[[inhabitant]] <- gene_tree[[inhabitant]] + branch$mu*minimum_time}
        else {gene_tree[[inhabitant]] <- branch$mu*minimum_time}
      }
      pair <- sample(branch$inhabitants, 2)
      new <- paste(pair[1], pair[2], sep = "")
      branch$inhabitants <- setdiff(c(branch$inhabitants, new), pair)
      inner <- solve_branch(branch, gene_tree, minimum_time + start, end)
      branch <- inner$branch
      gene_tree <- inner$gene_tree
    }
    else { # If there's no coalescence events, we only record the time passing.
      for (inhabitant in branch$inhabitants) {
        if (inhabitant %in% names(gene_tree)) {gene_tree[[inhabitant]] <- gene_tree[[inhabitant]] + branch$mu*(end - start)}
        else {gene_tree[[inhabitant]] <- branch$mu*(end - start)}
      }
    }
  }
  return(list(branch = branch, gene_tree = gene_tree))
}

solve_tree <- function(branches) {
  gene_tree <- list()
  # Compile a vector of pivotal moments (individual sampling times, population divergence and admixture times).
  pivotal <- numeric(0)
  for (branch in branches) {
    for (sample in branch$samples) {pivotal <- c(pivotal, sample)}
    for (switch in branch$switches) {pivotal <- c(pivotal, switch[1])}
  }
  pivotal <- unique(pivotal)
  pivotal <- pivotal[order(pivotal)]
  # Add infinity as the last pivotal moment.
  pivotal[length(pivotal) + 1] <- Inf
  # For all pivotal moments except the infinity:
  for (t in seq(1, length(pivotal) - 1)) {
    now <- pivotal[t]
    soon <- pivotal[t + 1]
    # First perform the special actions of the pivotal moment (samplings, switches).
    for (branch in names(branches)) {
      for (sample in names(branches[[branch]]$samples)) {
        if (branches[[branch]]$samples[[sample]] == now) {
          branches[[branch]]$inhabitants <- c(branches[[branch]]$inhabitants, sample)
        }
      }
      for (switch in names(branches[[branch]]$switches)) {
        if (branches[[branch]]$switches[[switch]][1] == now) {
          # Moving inhabitants from a branch to another, possibly.
          for (inhabitant in branches[[branch]]$inhabitants) {
            if (runif(1) < branches[[branch]]$switches[[switch]][2]) {
              branches[[branch]]$inhabitants <- setdiff(branches[[branch]]$inhabitants, inhabitant)
              branches[[switch]]$inhabitants <- c(branches[[switch]]$inhabitants, inhabitant)
            }
          }
        }
      }
    }
    # Then solve each branch until the next pivotal moment.
    for (branch in names(branches)) {
      solved <- solve_branch(branches[[branch]], gene_tree, now, soon)
      branches[[branch]] <- solved$branch
      gene_tree <- solved$gene_tree
    }
  }
  return(gene_tree)
}

pattern <- function(tree, Tmax = 0) {
  # Adding a number of point mutations somewhere in the gene tree and reporting the ABBABABA-pattern.
  L <- numeric(0)
  for (name in names(tree)) {
    if (tree[[name]] < Inf) {L[length(L) + 1] <- tree[[name]]}
    else {L[length(L) + 1] <- 0}
  }
  # We condition on a pattern being polymorphic, but often there's no mutation at all.
  # We use a resampling scheme that takes to total gene tree length and the possibility of multiple mutations into account.
  # The resampling can be made a bit faster as long as the total gene tree length 'sum(L)' can't exceed a predetermined maximum value 'Tmax'.
  if (Tmax == 0) {
    amount <- 1 # If we don't want to bother with bias correction and recurrent mutations we can just set 'Tmax' to 0 and assume a single mutation.
  }
  else {
    if (sum(L) > Tmax) {stop("The argument 'Tmax' must be chosen to be so big that the overall length of the gene tree virtually never exceeds it.")}
    amount_done <- FALSE
    while (amount_done == FALSE) {
      amount <- rpois(1, sum(L))
      if (amount == 0) {
        if (runif(1) > dpois(0, Tmax)/dpois(0, sum(L))) {amount_done <- TRUE} # In this case we need to randomize a new gene tree.
        # Otherwise we can randomize a new amount within the same gene tree.
      } else {amount_done <- TRUE} # A nonzero 'amount' is accepted.
    }
    if (amount == 0) {return(NA)}
  }
  mutations <- runif(amount)*sum(L)
  changes <- character(0) # The first nodes that are affected by the mutations.
  for (mutation in mutations) {changes[length(changes) + 1] <- names(tree)[sum(mutation > cumsum(L)) + 1]}
  # This will only work if populations are coded 'P1', 'P2', ... and there's less than ten of them.
  changes <- paste(changes, collapse = "")
  leaves <- max(nchar(names(tree)))/2 # Divided by two because of the character 'P'.
  pattern <- rep("A", leaves)
  for (j in seq(1, leaves)) {if (length(regmatches(changes, gregexpr(paste(j), changes))[[1]]) %% 2) {pattern[j] <- "B"}}
  if (length(unique(pattern)) > 1) {
    pattern <- paste(pattern, collapse = "")
    return(pattern)
  } # Returning a polymorphic pattern.
  else {return(NA)} # There was mutations but the resulting pattern was monomorphic.
}



## ----existing statistics------------------------------------------------------

n <- function(X, x) {
	if (x %in% names(X)) {return(X[x])} else {return(0)}
}

DPA <- function(X, alpha = 0.0027) { # The classic D-statistic.
	L <- n(X, "BABA") + n(X, "ABAB")
	R <- n(X, "ABBA") + n(X, "BAAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DFO <- function(X, alpha = 0.01) { # D_FO of Pease & Hahn.
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "BBBAA") + n(X, "AAABB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "AAABA") + n(X, "BBBAB")
	R <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "BBABA") + n(X, "AABAB") + n(X, "ABBAA") + n(X, "BAABB") + n(X, "AABAA") + n(X, "BBABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DIL <- function(X, alpha = 0.01) { # D_IL of Pease & Hahn.
	L <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "BBBAA") + n(X, "AAABB") + n(X, "BAABA") + n(X, "ABBAB") + n(X, "AAABA") + n(X, "BBBAB")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "BBABA") + n(X, "AABAB") + n(X, "BABAA") + n(X, "ABABB") + n(X, "AABAA") + n(X, "BBABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DFI <- function(X, alpha = 0.01) { # D_FI of Pease & Hahn.
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "BABBA") + n(X, "ABAAB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "ABAAA") + n(X, "BABBB")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "ABBBA") + n(X, "BAAAB") + n(X, "BAABA") + n(X, "ABBAB") + n(X, "BAAAA") + n(X, "ABBBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DOL <- function(X, alpha = 0.01) { # D_OL of Pease and Hahn.
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "BABBA") + n(X, "ABAAB") + n(X, "ABBAA") + n(X, "BAABB") + n(X, "ABAAA") + n(X, "BABBB")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "ABBBA") + n(X, "BAAAB") + n(X, "BABAA") + n(X, "ABABB") + n(X, "BAAAA") + n(X, "ABBBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

D1 <- function(X, alpha = 0.00333) { # D_1 of Eaton & Ree.
	L <- n(X, "ABBAA") + n(X, "BAABB")
	R <- n(X, "BABAA") + n(X, "ABABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

D2 <- function(X, alpha = 0.00333) { # D_2 of Eaton & Ree.
	L <- n(X, "ABABA") + n(X, "BABAB")
	R <- n(X, "BAABA") + n(X, "ABBAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

D12 <- function(X, alpha = 0.00333) { # D_12 of Eaton & Ree.
	L <- n(X, "BAAAB") + n(X, "ABBBA")
	R <- n(X, "ABAAB") + n(X, "BABBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}



## ----delta-statistics---------------------------------------------------------

n <- function(X, x) {
	if (x %in% names(X)) {return(X[x])} else {return(0)}
}

DS1 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB")
	R <- n(X, "BAABA") + n(X, "ABBAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS2 <- function(X, alpha = 0.01) {
	L <- n(X, "ABBAA") + n(X, "BAABB")
	R <- n(X, "ABABA") + n(X, "BABAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS3 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB")
	R <- n(X, "ABBAA") + n(X, "BAABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS4 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB")
	R <- n(X, "ABABA") + n(X, "BABAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS5 <- function(X, alpha = 0.01) {
	L <- n(X, "ABBBA") + n(X, "BAAAB")
	R <- n(X, "BABBA") + n(X, "ABAAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS6 <- function(X, alpha = 0.01) {
	L <- n(X, "BBABA") + n(X, "AABAB")
	R <- n(X, "BBBAA") + n(X, "AAABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS7 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS8 <- function(X, alpha = 0.01) {
	L <- n(X, "AABAA") + n(X, "BBABB")
	R <- n(X, "AAABA") + n(X, "BBBAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DA1 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB")
	R <- n(X, "ABBAA") + n(X, "BAABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DA2 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB")
	R <- n(X, "ABABA") + n(X, "BABAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DA3 <- function(X, alpha = 0.01) {
	L <- n(X, "ABBBA") + n(X, "BAAAB")
	R <- n(X, "BABBA") + n(X, "ABAAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DA4 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ1 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB")
	R <- n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ2 <- function(X, alpha = 0.01) {
	L <- n(X, "ABABA") + n(X, "BABAB")
	R <- n(X, "ABAAB") + n(X, "BABBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ3 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB")
	R <- n(X, "ABABA") + n(X, "BABAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ4 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAB") + n(X, "ABBBA")
	R <- n(X, "ABAAB") + n(X, "BABBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ5 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB")
	R <- n(X, "ABBAA") + n(X, "BAABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ6 <- function(X, alpha = 0.01) {
	L <- n(X, "AABBA") + n(X, "BBAAB")
	R <- n(X, "AABAB") + n(X, "BBABA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ7 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ8 <- function(X, alpha = 0.01) {
	L <- n(X, "AAABA") + n(X, "BBBAB")
	R <- n(X, "AAAAB") + n(X, "BBBBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS16 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "AAABB") + n(X, "BBBAA")
	R <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "AABAB") + n(X, "BBABA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS26 <- function(X, alpha = 0.01) {
	L <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "AAABB") + n(X, "BBBAA")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "AABAB") + n(X, "BBABA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS35 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "ABAAB") + n(X, "BABBA")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS45 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BABBA")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS57 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAB") + n(X, "ABBBA") + n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABAAB") + n(X, "BABBA") + n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS68 <- function(X, alpha = 0.01) {
	L <- n(X, "AABAB") + n(X, "BBABA") + n(X, "AABAA") + n(X, "BBABB")
	R <- n(X, "AAABB") + n(X, "BBBAA") + n(X, "AAABA") + n(X, "BBBAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS3457 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BABBA") + n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABBBA") + n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DS1268 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "ABBAA") + n(X, "BAABB") + n(X, "AAABB") + n(X, "BBBAA") + n(X, "AABAA") + n(X, "BBABB")
	R <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "AABAB") + n(X, "BBABA") + n(X, "AAABA") + n(X, "BBBAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DA12 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "ABABA") + n(X, "BABAB")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "BAABA") + n(X, "ABBAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DA13 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "ABAAB") + n(X, "BABBA")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DA23 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BABBA")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DA1234 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BABBA") + n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABBBA") + n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ16 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "AABAB") + n(X, "BBABA")
	R <- n(X, "BAAAB") + n(X, "ABBBA") + n(X, "AABBA") + n(X, "BBAAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ26 <- function(X, alpha = 0.01) {
	L <- n(X, "ABABA") + n(X, "BABAB") + n(X, "AABAB") + n(X, "BBABA")
	R <- n(X, "ABAAB") + n(X, "BABBA") + n(X, "AABBA") + n(X, "BBAAB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ35 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABBAA") + n(X, "BAABB")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "BABAA") + n(X, "ABABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ45 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAB") + n(X, "ABBBA") + n(X, "ABBAA") + n(X, "BAABB")
	R <- n(X, "ABAAB") + n(X, "BABBA") + n(X, "BABAA") + n(X, "ABABB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ57 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ68 <- function(X, alpha = 0.01) {
	L <- n(X, "AABBA") + n(X, "BBAAB") + n(X, "AAABA") + n(X, "BBBAB")
	R <- n(X, "AABAB") + n(X, "BBABA") + n(X, "AAAAB") + n(X, "BBBBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ3457 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "BAAAB") + n(X, "ABBBA") + n(X, "ABBAA") + n(X, "BAABB") + n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "ABAAB") + n(X, "BABBA") + n(X, "BABAA") + n(X, "ABABB") + n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}

DQ1268 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "AABAB") + n(X, "BBABA") + n(X, "AAABA") + n(X, "BBBAB")
	R <- n(X, "BAAAB") + n(X, "ABBBA") + n(X, "ABAAB") + n(X, "BABBA") + n(X, "AABBA") + n(X, "BBAAB") + n(X, "AAAAB") + n(X, "BBBBA")
	value <- as.numeric((L - R)/(L + R))
	zvalue <- (L - R)/sqrt(L + R)
	ci <- qnorm(1 - 0.5*alpha)*sqrt(4*L*R/(L + R)**3)
	pvalue <- 2*as.numeric(1 - pnorm(abs(zvalue)))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, zvalue, sign, ci))
}


