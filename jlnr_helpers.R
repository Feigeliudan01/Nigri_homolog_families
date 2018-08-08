##### HELP PROGRAM #####

################################################################################
# Install and load packages
################################################################################
# Remember to source the script before using the functions
source("https://bioconductor.org/biocLite.R")

# OBS: It is important to load the package: "cowplot" to keep the x,y axis
list.of.packages <- c("RMySQL", "data.table", "ape", "reshape", "ggplot2", "RColorBrewer",
                      "dplyr", "cluster", "geiger", "gtools", "scales", "gdata", "gplots",
                      "grid","gridExtra", "cowplot", "phylobase", "adephylo")

# Find packages not already installed - else install
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {biocLite(new.packages)}

# Load the packages
loads <- lapply(list.of.packages, library, character.only = TRUE)

# Install from source before use
#install.packages("devtools")
#devtools::install_github("lionel-/ggstance")
#devtools::install_github("GuangchuangYu/treeio")  
#devtools::install_github("GuangchuangYu/ggtree")  

library(devtools)
library(ggstance)
library(ggtree)

################################################################################
# PLOT CHANGES - BAR, PIE, PHYLOGENY
################################################################################

#---------------------------------------------------------#
# Blank background theme for ggplot2 - pie and bar charts
#---------------------------------------------------------#
# add to ggplot as a "+ black_theme"
blank_theme <- theme_minimal()+
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=14, face="bold")
    )


#---------------------------------------------------------#
# Costumized DTU colours
#---------------------------------------------------------#
total_colors <- function ()  {
    cbPalette <- c("#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#009E73")
    dtu_yellow_Palette <- c("#FFCC00", "#FF9900", "#FF6600", "#FF0000")
    dtu_red_Palette <- c("#990000", "#FF0099", "#CC3399", "#990066")
    dtu_blue_Palette <- c("#660066", "#660099", "#3366CC", "#33CCFF")
    dtu_green_Palette <- c("#99CC33", "#66CC00")
    dtu_prim_Palette <- c("#990000", "#999999")
    
    total_Pallet <- c(dtu_yellow_Palette, dtu_red_Palette, dtu_blue_Palette, dtu_green_Palette, dtu_prim_Palette)
    return(total_Pallet)
}

#---------------------------------------------------------#
# Rotate nodes in ape or nexus formatted tree
#---------------------------------------------------------#
# inputs: (tree , c(36,7,8) or "all")
rotateNodes<-function(tree,nodes,polytom=c(1,2),...){
    n<-length(tree$tip.label)
    if(nodes[1]=="all") nodes<-1:tree$Nnode+n
    for(i in 1:length(nodes)) 
        tree<-rotate(tree,nodes[i],polytom)
    if(hasArg(reversible)) reversible<-list(...)$reversible
    else reversible<-TRUE
    if(reversible){ 
        ii<-which(tree$edge[,2]<=n)
        jj<-tree$edge[ii,2]
        tree$edge[ii,2]<-1:n
        tree$tip.label<-tree$tip.label[jj]
    }
    return(tree)
}




################################################################################
# MYSQL CONNECTION AND DATA RETRIVAL
################################################################################

#---------------------------------------------------------#
# Provide query for mysql and get your result as dataframe, closes connection afterwards
#---------------------------------------------------------#
# input: query <- ("SELECT ..."), user_name, pasword, database_name, database_host. query is enough.
aspDbFetch <- function(query = NULL, user="remotejlnr_87", password="YY_4wazmbR", dbname="aspminedb", host="192.38.13.24", port = 3306){
    if(is.null(query)==FALSE){
        con <- dbConnect(MySQL(), user=user, password=password, dbname=dbname, host=host, port = port)
        #print(paste0("Connection established to ", dbname))
        aspData <- dbGetQuery(con,query)
        #print(sprintf("Fetched result with %s rows and %s columns", dim(aspData)[1],dim(aspData)[2]))
        dbDisconnect(con)
        #print("Closing connection")
        return(aspData)

    }else{
        print("Please check your query")
    }
}
