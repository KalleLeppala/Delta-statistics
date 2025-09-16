########## COMPUTING DELTA-STATISTICS FROM PLINK TRAW-FILES #########################################################

import math

# mymean
# Computes half (because diploid) the mean of the provided list.
#
# list	A list of strings that are allele counts of diploid individuals from a population or missing data coded as "NA".
#
# Returns the allele frequency unless all data was missing, in which case returns "NA".
def mymean(list):
    sum = 0
    amo = 0
    for item in list:
        if item != "NA":
            sum = sum + int(item)
            amo = amo + 1
    if amo == 0:
        return("NA")
    else:
        return(0.5*sum/amo)

# Delta
# Computes estimators of Delta-statistics and variances of those estimators.
#
# traw		A PLINK traw-file created with the command --recode A-transpose. Contains a header row and one row for each genetic locus.
#		First six columns are chromosome, ID, position in cmorgans, base-pair coordinate, counted allele and the reference allele.
#		Of these only the chromosome and the base-pair coordinate is used by the function.
#		The rest of the columns are allele counts of diploid individuals, missing data coded as "NA".
# populations	List of K lists. The K:th contain indices of the columns corresponding to the K:th population.
#		Keep in mind Python starts counting from zero and the first six columns contain other stuff, so the smallest possible index is 6.
# left		List of patterns that contribute positively to the Delta-statistic, written as words in {A, B}^K where B is the counted allele.
# right		List of patterns that contribute negatively to the Delta-statistic, written as words in {A, B}^K where B is the counted allele.
# blocksize	The size of blocks used, defaults at two millions.
#		A new block starts when the chromosome changes or the blocksize is reached in base pair coordinates.
#               To suppress the use of blocks (only when the parameter correction is False), set the blocksize to 0.
# correction	Depricated parameter that does nothing, defaults at False.
#               Program performs a window-wise computation of the statistics where no correction for LD is attempted, in the style of [Pease & Hahn 2015].
#		The statistic within a window is simply (L - R)/(L + R) and the Z-score is simply (L - R)/sqrt(L + R).
#               The program also prints on the console a single consensus statistic where inference is based in block jackknifing in the style of [Green 2010].
#		This block jackknifing  technique [Busin 1999] is described in a file called "wjack.pdf" by Nick Patterson you can find by Googling.
# output	The name of the file where results are printed, using one row per window if the parameter correction is False.
#		The default value is "" - then the results are not printed.
# progress	Whether the fuction prints progress to console or not. Defaults at True because I get restless otherwise.
#
# The function returns a list of lists of two numbers.
# The first one is the estimator of the Delta-statistic, the second is the Z-score i.e. the estimator divided by its estimated standard deviation.
# The program also prints on the console a consensus statistic and its Z-score that is based on variance estimated using a block jackknife.
def Delta(traw,
          populations,
          left,
          right,
          blocksize = 2000000,
	  correction = False,
	  output = "",
	  progress = True):
    f = open(traw, "r")
    f.readline() # Removing the header.
    counts = [] # This will contain one list [L, R, N] per block.
    chromosome = -1 # Impossible value.
    L = 0
    R = 0
    N = 0
    counter = 0
    for line in f:
        if counter % 1000 == 0 and progress == True:
            print("Lines read: " + str(counter))
        counter = counter + 1
        if blocksize != 0: # If we're not using blocks at all, we just skip the next few lines.
            # If a block is ready, we save results and start a new block.
            if int(line.split()[0]) != chromosome:
                chromosome = int(line.split()[0])
                position = int(line.split()[3])
                counts.append([L, R, N])
                L = 0
                R = 0
                N = 0
            elif int(line.split()[3]) > position + blocksize:
                position = int(line.split()[3])
                counts.append([L, R, N])
                L = 0
                R = 0
                N = 0
        # Treat the variable.
        missing = False
        for word in left:
            product = 1
            for i in range(len(populations)):
                frequency = mymean([line.split()[j] for j in populations[i]])
                if frequency == "NA":
                    product = 0
                    missing = True
                else:
                    if word[i] == "B":
                        product = product*frequency
                    else:
                        product = product*(1 - frequency)
            L = L + product
        for word in right:
            product = 1
            for i in range(len(populations)):
                frequency = mymean([line.split()[j] for j in populations[i]])
                if frequency == "NA":
                    product = 0
                    missing = True
                else:
                    if word[i] == "B":
                        product = product*frequency
                    else:
                        product = product*(1 - frequency)
            R = R + product
        if missing == False:
            N = N + 1
    f.close()
    counts.append([L, R, N])
    if blocksize != 0:
        counts = counts[1:] # The first element is just zeroes so we remove it.
    # Compute the estimators.
    result = []
    for i in range(len(counts)):
        if counts[i][0] + counts[i][1] > 0:
            result.append([(counts[i][0] - counts[i][1])/(counts[i][0] + counts[i][1]), (counts[i][0] - counts[i][1])/math.sqrt(counts[i][0] + counts[i][1])])
        else:
            result.append(["NA", "NA"])
    if output != "":
        g = open(output, "w")
        for window in result:
            nothing = g.write(str(window[0]) + " " + str(window[1]) + "\n")
    # Compute a consensus statistic using a block jackknife:
    D = [i[0] for i in result]
    D = [float(i) for i in D if i != "NA"]
    meanD = sum(D)/len(D)
    seD = 0
    for i in range(0, len(D)):
        sub = D[:i] + D[(i + 1):]
        seD = seD + (meanD - sum(sub)/len(sub))**2
    seD = math.sqrt((len(D) - 1)*seD/len(D))
    zD = meanD/seD
    print("Consensus statistic and its Z-score: " + str([meanD, zD]))
    return(result)

########## SHORTCUTS ################################################################################################

# D
# The classic Patterson's D-statistic [Green 2010].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def D(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABA", "ABAB"]
    right = ["ABBA", "BAAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DFO
# The D_FO-statistic from the D_FOIL -method [Pease 2015].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DFO(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "BBBAA", "ABABA", "AAABA", "ABABB", "AAABB", "BABAB", "BBBAB"]
    right = ["BAABA", "BBABA", "ABBAA", "AABAA", "ABBAB", "AABAB", "BAABB", "BBABB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DIL
# The D_IL-statistic from the D_FOIL -method [Pease 2015].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DIL(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["ABBAA", "BBBAA", "BAABA", "AAABA", "BAABB", "AAABB", "ABBAB", "BBBAB"]
    right = ["ABABA", "BBABA", "BABAA", "AABAA", "BABAB", "AABAB", "ABABB", "BBABB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DFI
# The D_FI-statistic from the D_FOIL -method [Pease 2015].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DFI(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "BABBA", "ABABA", "ABAAA", "ABABB", "ABAAB", "BABAB", "BABBB"]
    right = ["ABBAA", "ABBBA", "BAABA", "BAAAA", "BAABB", "BAAAB", "ABBAB", "ABBBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DOL
# The D_OL-statistic from the D_FOIL -method [Pease 2015].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAABA", "BABBA", "ABBAA", "ABAAA", "ABBAB", "ABAAB", "BAABB", "BABBB"]
    right = ["ABABA", "ABBBA", "BABAA", "BAAAA", "BABAB", "BAAAB", "ABABB", "ABBBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS16
# The ''freshman sum'' of DS1 and minus DS6.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DS16(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "AAABB", "BBBAA"]
    right = ["BAABA", "ABBAB", "AABAB", "BBABA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS26
# The ''freshman sum'' of DS2 and minus DS6.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DS26(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["ABBAA", "BAABB", "AAABB", "BBBAA"]
    right = ["ABABA", "BABAB", "AABAB", "BBABA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS35
# The ''freshman sum'' of DS3 and minus DS5.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DS35(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "ABAAB", "BABBA"]
    right = ["ABBAA", "BAABB", "BAAAB", "ABBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS45
# The ''freshman sum'' of DS4 and minus DS5.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DS45(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAABA", "ABBAB", "ABAAB", "BABBA"]
    right = ["ABABA", "BABAB", "BAAAB", "ABBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS57
# The ''freshman sum'' of DS5 and DS7.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DS57(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAAAB", "ABBBA", "BAAAA", "ABBBB"]
    right = ["ABAAB", "BABBA", "ABAAA", "BABBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS68
# The ''freshman sum'' of DS6 and DS8.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DS68(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["AABAB", "BBABA", "AABAA", "BBABB"]
    right = ["AAABB", "BBBAA", "AAABA", "BBBAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS3457
# The ''freshman sum'' of DS3, DS4, minus DS5 and DS7.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DS3457(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "BAABA", "ABBAB", "ABAAB", "BABBA", "BAAAA", "ABBBB"]
    right = ["ABBAA", "BAABB", "ABABA", "BABAB", "BAAAB", "ABBBA", "ABAAA", "BABBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS1268
# The ''freshman sum'' of DS1, DS2, minus DS6 and DS8.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DS1268(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "ABBAA", "BAABB", "AAABB", "BBBAA", "AABAA", "BBABB"]
    right = ["BAABA", "ABBAB", "ABABA", "BABAB", "AABAB", "BBABA", "AAABA", "BBBAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DA12
# The ''freshman sum'' of DA1 and minus DA2.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DA12(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "ABABA", "BABAB"]
    right = ["ABBAA", "BAABB", "BAABA", "ABBAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DA13
# The ''freshman sum'' of DA1 and minus DA3.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DA13(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "ABAAB", "BABBA"]
    right = ["ABBAA", "BAABB", "BAAAB", "ABBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DA23
# The ''freshman sum'' of DA2 and minus DA3.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DA23(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAABA", "ABBAB", "ABAAB", "BABBA"]
    right = ["ABABA", "BABAB", "BAAAB", "ABBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DA1234
# The ''freshman sum'' of DA1, DA2, minus DA3 and DA4.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DA1234(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "BAABA", "ABBAB", "ABAAB", "BABBA", "BAAAA", "ABBBB"]
    right = ["ABBAA", "BAABB", "ABABA", "BABAB", "BAAAB", "ABBBA", "ABAAA", "BABBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DQ16
# The ''freshman sum'' of DQ1 and minus DQ6.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DQ16(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAABA", "ABBAB", "AABAB", "BBABA"]
    right = ["BAAAB", "ABBBA", "AABBA", "BBAAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DQ26
# The ''freshman sum'' of DQ2 and minus DQ6.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DQ26(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["ABABA", "BABAB", "AABAB", "BBABA"]
    right = ["ABAAB", "BABBA", "AABBA", "BBAAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DQ35
# The ''freshman sum'' of DQ3 and minus DQ5.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DQ35(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAABA", "ABBAB", "ABBAA", "BAABB"]
    right = ["ABABA", "BABAB", "BABAA", "ABABB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DQ45
# The ''freshman sum'' of DQ4 and minus DQ5.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DQ45(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAAAB", "ABBBA", "ABBAA", "BAABB"]
    right = ["ABAAB", "BABBA", "BABAA", "ABABB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DQ57
# The ''freshman sum'' of DQ5 and DQ7.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DQ57(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "BAAAA", "ABBBB"]
    right = ["ABBAA", "BAABB", "ABAAA", "BABBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DQ68
# The ''freshman sum'' of DQ6 and DQ8.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DQ68(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["AABBA", "BBAAB", "AAABA", "BBBAB"]
    right = ["AABAB", "BBABA", "AAAAB", "BBBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DQ3457
# The ''freshman sum'' of DQ3, DQ4, minus DQ5 and DQ7.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DQ3457(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAABA", "ABBAB", "BAAAB", "ABBBA", "ABBAA", "BAABB", "BAAAA", "ABBBB"]
    right = ["ABABA", "BABAB", "ABAAB", "BABBA", "BABAA", "ABABB", "ABAAA", "BABBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DQ1268
# The ''freshman sum'' of DQ1, DQ2, minus DQ6 and DQ8.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DQ1268(traw, populations, blocksize = 2000000, correction = False, output = "", progress = True):
    left = ["BAABA", "ABBAB", "ABABA", "BABAB", "AABAB", "BBABA", "AAABA", "BBBAB"]
    right = ["BAAAB", "ABBBA", "ABAAB", "BABBA", "AABBA", "BBAAB", "AAAAB", "BBBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))