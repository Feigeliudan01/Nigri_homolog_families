# -------------------------------------------------------------------------
#                       ADDITIONAL IMPORTS
# -------------------------------------------------------------------------
# Install from source before use
#install.packages("devtools")
#devtools::install_github("lionel-/ggstance")
#devtools::install_github("GuangchuangYu/treeio")  
#devtools::install_github("GuangchuangYu/ggtree")


# -------------------------------------------------------------------------
#                            SETTINGS
# -------------------------------------------------------------------------
# Input programs, files and tables
source('../jlnr_helpers.R')
source('../adjToPhylo.R')
source('../NigriPaper_CPU_IPRcoverage_tree.R')

treeFile = '../NigriPaper_phylotree.nex' 
IPR_homoTable = "P17_Nigri_hfam_IPR"

# Main directory
mainDir = '../Genetic_diversity'
outputDir = 'PanCoreUnique_KeyNumbers'
setwd(file.path(mainDir, outputDir))

# Flags and values
createIPRtables = TRUE
write_tables = TRUE
tabletitleAppend = 'P17_R'
appendTo = "F"
valueSeparation = "comma"
plotCol = c("#8CD1C7FF", "#CCCCCC","#3e3e3fff", "#CCCCCC", "#F5B561FF", "#CCCCCC")  



# -------------------------------------------------------------------------
#                        INDIVIDUAL ANALYSIS
# -------------------------------------------------------------------------

# 38_species
fileName = '38_species'
org = c("Aspacu1", "Aspbru1", "Aspell1", "Aspeuc1", "Asphet1", "Asphom1", 
        "Asplac1", "Aspneo1", "Asppip1", "Aspsac1", "Aspscle1", "Aspscl1", 
        "Aspuva1", "Aspvad1", "Aspvio1", "Aspfij1", "Aspibe1", "Aspjap1", 
        "Aspcos1", "Aspfo1", "Aspac1", "Aspbr1", "Aspca3", "Aspfl1", "Aspfu1", 
        "Aspnid1", "Aspni7", "Aspni_DSM_1", "Aspni_NRRL3_1", "Aspni_bvT_1", 
        "Aspor1", "Aspph1", "Asptu1", "Aspka1_1", "Pench1", "Neucr2", "Aspwel1", 
        "Aspind1")

PanCoreUnique(IPR_homoTable = IPR_homoTable,  org = org, write_tables = write_tables, 
                          valueSeparation = valueSeparation, fileName = fileName, outputDir = outputDir, 
                          mainDir = mainDir, treeFile = treeFile, createIPRtables = createIPRtables, 
                          tabletitleAppend = tabletitleAppend, appendTo = appendTo, 
                          plotCol = plotCol)

# 36_Aspergilli
fileName = '36_Aspergilli'
org = c("Aspacu1", "Aspbru1", "Aspell1", "Aspeuc1", "Asphet1", "Asphom1", 
        "Asplac1", "Aspneo1", "Asppip1", "Aspsac1", "Aspscle1", "Aspscl1", 
        "Aspuva1", "Aspvad1", "Aspvio1", "Aspfij1", "Aspibe1", "Aspjap1", 
        "Aspcos1", "Aspfo1", "Aspac1", "Aspbr1", "Aspca3", "Aspfl1", "Aspfu1", 
        "Aspnid1", "Aspni7", "Aspni_DSM_1", "Aspni_NRRL3_1", "Aspni_bvT_1", 
        "Aspor1", "Aspph1", "Asptu1", "Aspka1_1", "Aspwel1", "Aspind1")
length(org)

PanCoreUnique(IPR_homoTable = IPR_homoTable,  org = org, write_tables = write_tables, 
              valueSeparation = valueSeparation, fileName = fileName, outputDir = outputDir, 
              mainDir = mainDir, treeFile = treeFile, createIPRtables = createIPRtables, 
              tabletitleAppend = tabletitleAppend, appendTo = appendTo, 
              plotCol = plotCol)


# Nigri_section
fileName = 'Nigri_section'
org = c("Aspacu1", "Aspbru1", "Aspell1", "Aspeuc1", "Asphet1", "Asphom1", 
        "Asplac1", "Aspneo1", "Asppip1", "Aspsac1", "Aspscle1", "Aspscl1", 
        "Aspuva1", "Aspvad1", "Aspvio1", "Aspfij1", "Aspibe1", "Aspjap1", 
        "Aspcos1", "Aspfo1", "Aspac1", "Aspbr1", "Aspca3", "Aspni7", 
        "Aspni_DSM_1", "Aspni_NRRL3_1", "Aspni_bvT_1", "Aspph1", "Asptu1", 
        "Aspka1_1", "Aspwel1", "Aspind1")

PanCoreUnique(IPR_homoTable = IPR_homoTable,  org = org, write_tables = write_tables, 
              valueSeparation = valueSeparation, fileName = fileName, outputDir = outputDir, 
              mainDir = mainDir, treeFile = treeFile, createIPRtables = createIPRtables, 
              tabletitleAppend = tabletitleAppend, appendTo = appendTo, 
              plotCol = plotCol)

# Tubigensis_clade
fileName = 'Tubigensis_clade'
org = c("Aspneo1", "Asptu1", "Aspcos1", "Aspeuc1", "Aspfo1", "Aspka1_1", 
        "Asppip1", "Aspvad1")

PanCoreUnique(IPR_homoTable = IPR_homoTable,  org = org, write_tables = write_tables, 
              valueSeparation = valueSeparation, fileName = fileName, outputDir = outputDir, 
              mainDir = mainDir, treeFile = treeFile, createIPRtables = createIPRtables, 
              tabletitleAppend = tabletitleAppend, appendTo = appendTo, 
              plotCol = plotCol)

# Niger_clade
fileName = 'Niger_clade'
org = c("Aspni7", "Aspni_NRRL3_1", "Asplac1", "Aspni_bvT_1", "Aspni_DSM_1", 
        "Aspph1") # 

PanCoreUnique(IPR_homoTable = IPR_homoTable,  org = org, write_tables = write_tables, 
              valueSeparation = valueSeparation, fileName = fileName, outputDir = outputDir, 
              mainDir = mainDir, treeFile = treeFile, createIPRtables = createIPRtables, 
              tabletitleAppend = tabletitleAppend, appendTo = appendTo, 
              plotCol = plotCol)

fileName = 'Niger_clade_incl_Aspwel1'
org = c("Aspni7", "Aspni_NRRL3_1", "Asplac1", "Aspni_bvT_1", "Aspni_DSM_1", 
        "Aspph1", "Aspwel1") 

PanCoreUnique(IPR_homoTable = IPR_homoTable,  org = org, write_tables = write_tables, 
              valueSeparation = valueSeparation, fileName = fileName, outputDir = outputDir, 
              mainDir = mainDir, treeFile = treeFile, createIPRtables = createIPRtables, 
              tabletitleAppend = tabletitleAppend, appendTo = appendTo, 
              plotCol = plotCol)

