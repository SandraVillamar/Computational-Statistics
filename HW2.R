#PART 3
unexp <- c(8,11,12,14,20,43,111)
exp <- c(35,56,83,92,128,150,176,208)
expNew <- exp - 25
#HT
wilcox.test(unexp, expNew, alternative = "less", exact = TRUE)
#CI
wilcox.test(unexp, exp, exact = TRUE, conf.int = T)


#PART 4
pH <- c(7.02,7.34,7.28,7.09,7.45,7.40,7.32)
wilcox.test(pH, exact = TRUE, conf.int = T)

