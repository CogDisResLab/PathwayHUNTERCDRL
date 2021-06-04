#Stopwords
mystopwords <- tibble(word = c("frequency", "rate", "pathways","involving","chemical","reactions","assembly"))
#Data Analysis for Wilma
#Loading in data files from Wilma
mddatp <- read.csv("MDD ATPase extracted gene lists_p0.05_enrichment.csv")
mddnaka <- read.csv("MDD Nakamura sublibrary significant gene extracted list 0.05 enrichment.csv")
sczdeg <- read.csv("SCZ_DEG_A_p0.05_enrichment.csv")
shortnaka <- read.csv("short_Nakamura_ATP_SCZ_0.05_enrichment.csv")

#data clean up
#removing first 2 columns & renaming a the GOID column
mddatp <- subset(mddatp, select= -c(X,X.1))
names(mddatp)[1] <- "GOID"

mddnaka <- subset(mddnaka, select= -c(X,X.1))
names(mddnaka)[1] <- "GOID"

PathwayHUNTER(mddnaka, 8, "Test.csv", 2)
