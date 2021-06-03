#pathwayHUNTER Script
#Liraries needed for pathwayHUNTER
library(GeneOverlap)
library(dplyr)
library(GO.db)
library(git2r)
library(GO.db)
library(tibble)
library(vctrs)
library(tidyr)
library(stringr)
library(tm)
library(text2vec)
library(stopwords)
library(tidytext)
library(ggplot2)
library(ComplexHeatmap)


PathwayHUNTER <- function(x,y,z,num){
  
  #Loading in the big dataframe
  #This big dataframe contains all the cuts
  #We use variable y to cute the datafram into our 50 - 2250 clusters
  Bigt <- test[[y]]
  Bigt$level <- NULL
  my_list <- split(Bigt$goid,Bigt$cluster)
  
  #Grabbing the GOIDs from User input
  #Make it into a list so that we can use GeneOverlap package to compare aginst our Large dataframe
  Sgoids <- x$GOID
  sgoidlist <- list(c(Sgoids))
  gs.RNASeq <- 44509
  gom.obj <- newGOM(my_list,sgoidlist,gs.RNASeq)
  
  #Grabbing the pval of the overlap
  #nrow grabs the the number of rows
  #then we filter the clusters based on the .05 pval
  pvalmatrix <- base::as.data.frame(GeneOverlap::getMatrix(gom.obj, name="pval"))
  rows <- nrow(pvalmatrix)
  colnames(pvalmatrix) <- "pVal"
  pvalmatrix$Clusters <- 1:rows
  SubPval<- pvalmatrix %>% dplyr::filter(pVal< 1.0)
  
  #Based on the Clusters that survive our cut off we set that to xx
  #we then use pdata that contains all the clusters to GOID at each level cut off and set that to value
  xx <- SubPval$Clusters
  value <- pdata %>% dplyr::filter(level == rows, cluster %in%c(xx))
  UserValue <- merge(x,value,by.x = "GOID", by.y = "goid")
  UserValue <- as.data.frame(UserValue$GOID)
  names(UserValue)[1] <- "GOID"
  #This is the starts of the textmining/theme generating portion of the code
  #Merging GOID Clusters with Defeinitions
  #Grabbing GOID, TERM and Defeinition from GO.db package to merge with our GOID set
  keys <- (keys(GO.db::GO.db))
  GO.TERMS <- AnnotationDbi::select(GO.db,keys=keys, columns =c("GOID","TERM", "DEFINITION"))
  
  #Merging our GO.db package dataframe with our GOIDs that belong to individual clusters
  ClusterAnalysisWithDef <- merge(value, GO.TERMS, by.x = "goid", by.y= "GOID" )
  UserClusterAnalysisWithDef<- merge(UserValue, GO.TERMS, by.x = "GOID", by.y = "GOID")
  #This line replaces any NA values within the definition column with the Term column
  #This essential gives us better theme generating
  ClusterAnalysisWithDef$DEFINITION <- ifelse(is.na(ClusterAnalysisWithDef$DEFINITION), ClusterAnalysisWithDef$TERM, ClusterAnalysisWithDef$DEFINITION)
  #We now have if then statments that allow us to generate themes based on 2 word or one word analysis
  if (num ==2){ tidy_Clusters <- ClusterAnalysisWithDef %>%
    tidytext::unnest_tokens(word,DEFINITION, token = "ngrams", n= 2)
  Usertidy_Clusters <- UserClusterAnalysisWithDef %>%
    tidytext::unnest_tokens(word,DEFINITION, token = "ngrams", n= 2)
  
  #Separting the bi gram into 2 columns
  #Dividing the 2 part pharases into 2 sperate columns to be cleaned
  tidy_ClustersSep <- tidy_Clusters %>%
    tidyr::separate(word, c("word1", "word2"), sep= " ")
  Usertidy_ClustersSep <- Usertidy_Clusters %>%
    tidyr::separate(word, c("word1", "word2"), sep= " ")
  
  #Filtering the 2 separate columns by stopwords
  #Removal of common stop_words and our stopwords
  tidy_ClustersFilter <- tidy_ClustersSep %>%
    filter(!word1 %in% stop_words$word) %>%
    filter(!word2 %in% stop_words$word) %>%
    filter(!word1 %in% mystopwords$word) %>%
    filter(!word2 %in% mystopwords$word)
  Usertidy_ClustersFilter <- Usertidy_ClustersSep %>%
    filter(!word1 %in% stop_words$word) %>%
    filter(!word2 %in% stop_words$word) %>%
    filter(!word1 %in% mystopwords$word) %>%
    filter(!word2 %in% mystopwords$word)
  
  
  #unite the bigram
  tidycluster_unite <- tidy_ClustersFilter %>%
    tidyr::unite(word, word1, word2, sep = " ")
  Usertidycluster_unite <- Usertidy_ClustersFilter %>%
    tidyr::unite(word, word1, word2, sep = " ")}
  
  else{
      tidy_Clusters <- ClusterAnalysisWithDef %>%
        tidytext::unnest_tokens(word,DEFINITION, token = "ngrams", n= 1)
      
      #Separting the bi gram into 2 columns
      #Dividing the 2 part pharases into 2 sperate columns to be cleaned
      tidy_ClustersSep <- tidy_Clusters %>%
        tidyr::separate(word, c("word1"), sep= " ")
      
      #Filtering the 2 separate columns by stopwords
      #Removal of common stop_words and our stopwords
      tidy_ClustersFilter <- tidy_ClustersSep %>%
        filter(!word1 %in% stop_words$word) %>%
        filter(!word1 %in% mystopwords$word)
      
      
      #unite the bigram
      tidycluster_unite <- tidy_ClustersFilter %>%
        tidyr::unite(word, word1, sep = " ")}
  #Calculated tf, idf, and tf_idf by the bigram
  Themes<- tidycluster_unite %>%
    dplyr::count(cluster,word,sort=TRUE) %>%
    tidytext::bind_tf_idf(word,cluster,n) %>%
    dplyr::group_by(cluster) %>%
    top_n(1)
  #Themes2.0 is just where we have multiple / themes
  Themes2.0<- Themes %>% dplyr::group_by(cluster) %>% dplyr::mutate(ClusterTheme = paste0(word, collapse = " \\\ ")) %>%
    dplyr::distinct(cluster, .keep_all = T) %>% ungroup()
  
  Usertidycluster_unite <- merge(Usertidycluster_unite, value, by.x= "GOID", by.y= "goid")
  UserThemes <- Usertidycluster_unite %>%
    dplyr::count(cluster,word,sort=TRUE) %>%
    tidytext::bind_tf_idf(word,cluster,n) %>%
    dplyr::group_by(cluster) %>%
    top_n(1)
  
  UserThemes2.0<- UserThemes %>% dplyr::group_by(cluster) %>% dplyr::mutate(UserTheme = paste0(word, collapse = " \\\ ")) %>%
    dplyr::distinct(cluster, .keep_all = T) %>% ungroup()
  
  #Making the heatmap
  #Merge our original file with Enrichment scores with that of goid and cluster assignment
  #I then merge that with the themes generated for each cluste from our textmining
  myData <- base::merge(x, value, by.x = "GOID", by.y = "goid")
  myData <- base::merge(myData,Themes2.0, by.x ="cluster", by.y = "cluster")
  myData <- base::merge(myData,UserThemes2.0, by.x ="cluster", by.y = "cluster")
  myData[is.na(myData)] = 0
  
  #Marissa provided this
  p <- dplyr::select(myData, -c(level,n.x,tf.x,idf.x,tf_idf.x,n.y,tf.y,idf.y,tf_idf.y,word.x,word.y,GOID,ClusterWord,cluster))
  
  #Make a data matrix to work with Complex heatmap
  pp <- as.matrix(p)
  
  #Split to assign the left hand column of the heatmap the themes
  split <- paste0( myData$word)
  
  #This is still a work in progress but generating the heatmap
  #default.hmap <- ComplexHeatmap::Heatmap(pp,cluster_row_slices = FALSE,row_title_gp = gpar(font = .0001),
                                         # row_title_rot = 0,row_gap = unit( 4, "mm"),split = split,use_raster = F)
  
  #pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
  #ComplexHeatmap::draw(default.hmap, newpage=TRUE)
  
  #making table
  #Varibal z is how you can make a table just put the file name and it will make a table with
  Table <- dplyr::select(myData, -c(level,n.x,tf.x,idf.x,tf_idf.x,n.y,tf.y,idf.y,tf_idf.y,word.x,word.y,GOID,ClusterWord,cluster))
  Table <- merge(Table, SubPval,by.x = "cluster", by.y = "Clusters")
  # select variables v1, v2, v3
  myvars <- c("GOID", "TERM")
  TERMS <- GO.TERMS[myvars]
  Table <- merge(Table, TERMS, by.x = "GOID", by.y = "GOID")
  write.csv(Table, file= z)
}
