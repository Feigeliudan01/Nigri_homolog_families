#setwd('/home/seth/seth-1')

#source('seths_helpers.R')

adjToPhylogeny <- function(dfWithNames, treeFile, reverse_order = TRUE){
   print('ATTENTION: I changed this function so it takes real_name as input')
  # Input is a dataframe with a name column from our aspmine and a treefile
  # The treefile is used to sort the entries of the dataframe
  # I sort the whole dataframe and not only a vector because
  # I need this for my stacked barplots, please let me know if you need it for vectors as well
  
  #   > head(someDf)
  #       name  sm_short type
  #   1 Aspfo1      DMAT    2
  #   2 Aspfo1    HYBRID    6
  #   3 Aspfo1      NRPS   18
  #   4 Aspfo1 NRPS-Like   26
  #   5 Aspfo1       PKS   37
  #   6 Aspfo1  PKS-Like    6
  
  # The only thing you really need for the function to run is the name column
  
  library(dplyr)
  library(gdata)
  
  # Importing your treefile
  if(class(treeFile) == 'character'){
    fullTree <- read.tree(treeFile)
  } else if(class(treeFile) == 'phylo') {
    fullTree <- treeFile
    phyloNames <- fullTree$tip.label
  } else if(class(treeFile) %in% c('hclust', 'dendrogram')) {
    phyloNames <- labels(treeFile)
    } else {print('Please provide a valid input object of type character to load a file or type file to load the object directly')}
 
  # Cutting away names that are not included in our dfWithNames
  phyloNames <- intersect(phyloNames, dfWithNames$real_name)
  
  # Cutting away names that are not included in the df
  if(length(levels(dfWithNames$real_name)) > length(phyloNames)){
    warning('There are more names in the dataframe than in the phylogenetic tree. /n
          Reducing the dataframe names to names occuring in phylogenetic tree.')
    print('I will throw out:')
    print(setdiff(dfWithNames$real_name, phyloNames))
    dfWithNames <- dfWithNames[dfWithNames$real_name %in% phyloNames,]
  }
  
  phyloNames <- as.factor(phyloNames)
  dfWithNames$real_name <- as.factor(dfWithNames$real_name)

  if(reverse_order == TRUE){
    orderDummy <- rev(phyloNames)
  } else {orderDummy <- phyloNames}
  
  dfWithNames$real_name <- reorder.factor(dfWithNames$real_name, new.order = orderDummy)
  
  dfWithNames <- dfWithNames %>% arrange(real_name)
  
  return(dfWithNames)
  
}

# Trying out the function

#tree <-'RAxML_bestTree.uniq_qc'

#someDf<- aspDbFetch("select name, sm_short, count(*) as type from antismash join organism using (org_id) where sm_short != 'none' group by org_id, sm_short;")

#newDf <- adjToPhylogeny(dfWithNames = someDf, treeFile = tree)
