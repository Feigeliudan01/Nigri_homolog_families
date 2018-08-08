# Install from source before use
#install.packages("devtools")
#devtools::install_github("lionel-/ggstance")
#devtools::install_github("GuangchuangYu/treeio")  
#devtools::install_github("GuangchuangYu/ggtree")  

PanCoreUnique <- function(IPR_homoTable = NULL,  org = NULL, write_tables = FALSE, 
                          valueSeparation = "comma", fileName = 'test', outputDir = 'test', 
                          mainDir = "~/Desktop", treeFile = NULL, createIPRtables = FALSE, 
                          tabletitleAppend = '', appendTo = "front", 
                          plotCol = c("#249dc5ff", "#CCCCCC","#3e3e3fff", "#CCCCCC", "#e6a022ff", "#CCCCCC")){
    
    # -------------------------------------------------------------------------
    # CREATE AND SET DIRECTORIES
    # -------------------------------------------------------------------------
    oldwd <- getwd()
    dir.create(file.path(mainDir), showWarnings = FALSE)
    dir.create(file.path(mainDir, outputDir), showWarnings = FALSE)
    dir.create(file.path(mainDir, outputDir, fileName), showWarnings = FALSE)
    setwd(file.path(mainDir, outputDir, fileName))
    
    
    # -------------------------------------------------------------------------
    #   CHECK GIVEN INPUTS   
    # -------------------------------------------------------------------------
    # Check inputs are present
    if(is.null(org)==TRUE) {
        stop("# ERROR: Please specify either your organism(s).")}
    
    if(is.null(IPR_homoTable)==TRUE) {
        stop("# ERROR: Please specify the IPR homoPF table to use.")}
    
    if(is.null(treeFile)) {
        stop("# ERROR: Please specify the tree file to use (nexus format).")}
    
    if(fileName =="test") {
        print("# INFO: Your file names will be called 'test'")}
    
    if(outputDir =="test") {
        print("# INFO: Your output directory will be called 'test'")}
    
    if(mainDir =="~/Desktop") {
        print("# INFO: Your output directory will be saved on your desktop")}
    
    
    # Convert appendTo input to either "f" or "b"
    appendTo = tolower(appendTo)
    front_back = c("front", "f", "back", "b")
    
    if(!(appendTo %in% front_back)){
        stop(sprintf("# ERROR: You did not specify if you would like to append '%s' to the front or back of the table name.\n# Error: Please choose ('front' or 'back') for the 'appendTo' variable", tabletitleAppend))
    } else {
        if(appendTo == "front"){
            appendTo = "f"
        } else if(appendTo == "back"){
            appendTo = "b"
        }
    }
    
    # Convert appendTo input to either "f" or "b"
    valueSeparation = tolower(valueSeparation)
    dot_comma = c("comma", "coma", "c", "d", "dot")
    
    if(!(valueSeparation %in% dot_comma)){
        stop(sprintf("# ERROR: You did not specify if if the numbers should be comma or dot seperated.\n# ERROR: Please choose ('dot' or 'comma') for the 'valueSeparation' variable", tabletitleAppend))
    } else {
        if(valueSeparation == "comma"){
            valueSeparation = "c"
        } else if(valueSeparation == "coma"){
            valueSeparation = "c"
        } else if(valueSeparation == "dot"){
            valueSeparation = "d"
        }
    }
    
    # Check if tables exists
    tables <- c(IPR_homoTable)
    
    for (i in 1:length(tables)){
        tableCheck_query <- sprintf("Show tables LIKE '%s'", tables[i])
        tableCheck <- aspDbFetch(tableCheck_query)
        if (nrow(tableCheck) != 1) {
            stop(sprintf("\n# ERROR: The table %s was found %s times, please specify your table again.", tables[i], nrow(tableCheck) ))
        }
    }
    
    # Check if colors exists
    # Subfunction
    isColor <- function(x) {
        res <- try(col2rgb(x),silent=TRUE)
        return(!"try-error"%in%class(res))
    }
    
    
    if (length(plotCol) != 6) {
        stop(sprintf("# ERROR: You picked %s colors: '%s'. \n# ERROR: Please pick 4 different colours in this order: c('no_IPR_coverage', 'Pan', 'Core', 'Unique').", length(plotCol), paste(plotCol, collapse="', '")))
    } else {
        for (colour in plotCol) {
            if (identical(isColor(colour), FALSE)){
                stop(sprintf("\n# ERROR: The color '%s' does not exist. Please pick a new color.", colour))
            }
        }
    }
    
    
    ###########################################################################
    ################## SUB FUNCTIONS FOR BAP PLOTS ############################
    ###########################################################################
    # -------------------------------------------------------------------------
    # Create total barplots for pan/core/unique for proteins and protein families (hfams)
    # -------------------------------------------------------------------------
    prot_hfam_PanCoreUnique <- function(data, xcol, ycol, varcol, x_axis_label, valueSeparation, colorChoice) {
        max_value <- max(data$value) + 25000
        p <- ggplot(data, aes(x=xcol, y=ycol, fill = varcol, width=0.8))+ 
            geom_bar(stat="identity", colour="black", position = "dodge")+ 
            xlab(x_axis_label)
        
        # Converting numbers to include dot or comma 1000-separation
        if (valueSeparation == "d"){
            p <- p + scale_y_continuous(name = "", labels=function(x) format(x, decimal.mark=",", big.mark=".", scientific = FALSE), expand = c(0,0), limits = c(0,max_value)) #DOT
        } else {
            p <- p + scale_y_continuous(name = "", labels=function(x) format(x, decimal.mark=".", big.mark=",", scientific = FALSE), expand = c(0,0), limits = c(0,max_value)) #COMMA
        }
        
        p <- p + scale_fill_manual(name = '', labels = c('Pan genome  ', "Core core   ", 'Species unique genes'),
                                   values = colorChoice) +
            geom_text(aes(y=ycol, ymax=ycol, label=dot), 
                      position= position_dodge(width=0.8), vjust=-0.5, color="black", hjust = 0.5)+
            theme(axis.title.y = element_blank(),
                  axis.title.x = element_text(face="bold", colour="black", size=16),
                  axis.text.x  = element_text(colour="black",  vjust=0.5, size=14),
                  axis.text.y  = element_text(colour="black", vjust=0.5, size=14)) +
            theme(panel.background = element_blank(), axis.ticks = element_line(colour = "black"))  +
            theme(axis.line = element_line(colour = "black", size = 0.5)) +
            theme(legend.position="top") + 
            guides(fill = guide_legend(override.aes = list(colour = NULL), keywidth = 1, keyheight = 1)) +
            theme(legend.key = element_rect(colour = "black"))
        return(p)
    }
    
    # -------------------------------------------------------------------------
    # Create tree + Pan, core and unique bar plots incl. interPro coverage
    # -------------------------------------------------------------------------
    IPRcov_PanCoreUnique_perSpecies <- function(totalPANmelt, totalCOREmelt, totalUNIQUEmelt, tree, tree_width, plotCol) {
        
        # Calculate sum
        Ptotals <- totalPANmelt %>%
            group_by(Species) %>%
            summarize(total = sum(value))
        
        Ctotals <- totalCOREmelt %>%
            group_by(Species) %>%
            summarize(total = sum(value))
        
        Utotals <- totalUNIQUEmelt %>%
            group_by(Species) %>%
            summarize(total = sum(value))
        
        # Combine plots and data labels
        pan_plot <- facet_plot(tree, panel = 'Genome size', data = totalPANmelt,
                               geom = geom_barh,
                               mapping = aes(x = value, fill = as.factor(variable)),
                               stat='identity', color = "black", size = 0.25) 
        pan_label <- facet_plot(pan_plot, geom=geom_text, 
                                mapping=aes(x=total+max(total)*0.06, label=format(total, decimal.mark=".", big.mark=",", trim=TRUE)), 
                                data=Ptotals, panel='Genome size')
        core_plot <- facet_plot(pan_label, panel = 'Core genome', data = totalCOREmelt,
                                geom = geom_barh,
                                mapping = aes(x = value, fill = as.factor(variable)),
                                stat='identity', color = "black", size = 0.25)
        core_label <- facet_plot(core_plot, geom=geom_text, 
                                 mapping=aes(x=total+max(total)*0.06, label=format(total, decimal.mark=".", big.mark=",", trim=TRUE)), 
                                 data=Ctotals, panel='Core genome')
        unique_plot <- facet_plot(core_label, panel = 'Species unique genes', data = totalUNIQUEmelt,
                                  geom = geom_barh,
                                  mapping = aes(x = value, fill = as.factor(variable)),
                                  stat='identity', color = "black", size = 0.25)
        unique_label <- facet_plot(unique_plot, geom=geom_text, 
                                   mapping=aes(x=total+max(total)*0.05,label=format(total, decimal.mark=".", big.mark=",", trim=TRUE)), 
                                   data=Utotals, panel='Species unique genes') +
            scale_fill_manual(name = '', labels = c("Not annotated    ", 'InterPro annotated'), values = plotCol) +
            xlim_tree(tree_width) + #theme(strip.text = element_text(size = rel(2)))
            xlim_expand(max(Ctotals$total), 'Core genome') # Enlarge specific panels x-axis 
        
        return(unique_label)
    }
    
    # -------------------------------------------------------------------------
    # Create tree + Pan, core and unique bar plots
    # -------------------------------------------------------------------------
    PanCoreUnique_perSpecies <- function(PCU_melt, Ptotals, phy, plotCol, tree_width) {
        CUA_plot <- facet_plot(phy, panel = 'Core, unique and accessory seqments per genome', data = PCU_melt,
                               geom = geom_barh,
                               mapping = aes(x = value, fill = as.factor(variable)),
                               stat='identity', color = "black", size = 0.25) 
        CUA_label <- facet_plot(CUA_plot, geom=geom_text, mapping=aes(x=Pan+max(Pan)*0.06, 
                                                                      label=format(Pan, decimal.mark=".", big.mark=",", trim=TRUE)), 
                                data=Ptotals, panel='Core, unique and accessory seqments per genome') +
            scale_fill_manual(name = '', values = plotCol[c(3,5,1)]) + #c( "#a3a1b3ff", "#F5B561FF","#8CD1C7FF")
            xlim_tree(tree_width) + xlim_expand(max(Ptotals$Pan), 'Core, unique and accessory seqments per genome')
        
        return(CUA_label)
    }
    ###########################################################################
    ############################ MAIN PROGRAM #################################
    ###########################################################################
    
    ###########################################################################
    # RETRIEVE SPECIES NAMES FROM DATABASE 
    ###########################################################################
    
    # -------------------------------------------------------------------------
    # Fetch org_names from organism table
    # -------------------------------------------------------------------------
    # Retrieving the organims names that will be used for the analysis
    org_short <- paste(org, collapse="','") # convert input org_names
    org_querystring <- sprintf("SELECT name, concat(substring(genus, 1, 1), '. ', real_name)\
                               as real_name FROM organism WHERE name IN ('%s');", org_short)
    org_table <- aspDbFetch(org_querystring)
    
    # Check that all org_names are found in "organism" table 
    if(length(org) != nrow(org_table)) {
        stop(sprintf("# ERROR: %s of the organisms are not found in the MySQL organism table, please check your organisms for: %s .", length(org)-nrow(org_table),setdiff(org,org_table$name)))
    }
    
    # Converting org_names to vectors
    organisms <- as.vector(org_table$name)
    
    # -------------------------------------------------------------------------
    # Check presense of org_names in homoFP table
    # -------------------------------------------------------------------------
    # Collect all organism names from tables
    homoTable_orgs_query <- sprintf("SELECT DISTINCT org_name from %s;", IPR_homoTable)
    homoTable_orgs <- aspDbFetch(homoTable_orgs_query)[[1]]
    
    # Collect non-found organisms in tables and exit if they exists
    not_found_org_homo = character(0)
    for (i in 1:length(organisms)){
        if(organisms[i] %in% homoTable_orgs == FALSE) {
            not_found_org_homo  <- c(not_found_org_homo, organisms[i])
        }
        
        # Exit if non-found organisms exists
        if(length(not_found_org_homo) != 0){
            homo_rename <- paste(not_found_org_homo, collapse=", ")
            stop(sprintf("# ERROR: %s are not found in %s, please check your organisms and table.", homo_rename, IPR_homoTable))
        }
    }
    
    ###########################################################################
    # ADJUST TREE FILE 
    ###########################################################################
    # Read tree file and check if input orgs are in the tree file
    print("# INFO: It is important that the tree file is in the right order for the program to work.")
    fullTree <- read.tree(treeFile)
    
    if(all(organisms %in% fullTree$tip.label) != TRUE) {
        stop(sprintf("# ERROR: The tree file doesn't contain the organism [%s], please check your organism input.", organisms[!organisms %in% fullTree$tip.label]))
    }
    
    # Subset tree to only contain input organisms
    dropOuters <- fullTree$tip.label[!fullTree$tip.label %in% organisms]
    if(!identical(dropOuters, character(0))){
        print('# INFO: Modifying tree for better overview')
        fullTree <- drop.tip(fullTree, dropOuters)
    }
    
    # Retrieve the order
    # IMPORTANT FOR SORTING THE OUTPUT FILES
    # OBS: THIS DOES NOT WORK
    is_tip <- fullTree$edge[,2] <= length(fullTree$tip.label)
    ordered_tips <- fullTree$edge[is_tip, 2]
    name_ordered_tips <- fullTree$tip.label[ordered_tips]
    
    # Rename 
    fullTree_newNames <- fullTree
    newNames <- lapply(fullTree$tip.label, function(x) {
        org_table[match(x, org_table$name),2]
        
    })
    newNames <- as.character(newNames)
    fullTree_newNames$tip.label <- newNames
    
    
    ###########################################################################
    # CREATE AND INDEX TEMP_IPR_TABLE - Incl. only the species of interest
    ###########################################################################
    # Including all homolog families and unique paralog families with their proteins and interpro domains
    IPRtempTable <- sprintf("t_RtempIPR_%s", fileName)
    checkTable <- aspDbFetch(sprintf("SHOW TABLES LIKE '%s';", IPRtempTable))
    
    if(nrow(checkTable)==0){
        print('# INFO: Creating temporary InterPro table with the selected organisms')
        organism_formattet <- toString(shQuote(mixedsort(organisms)))
        IPRtempTable_query <- sprintf("CREATE TABLE %s
                                      SELECT a.hfam, org_name, protein_id, ipr_id, b.org_count
                                      FROM %s a
                                      JOIN
                                      (SELECT hfam, COUNT(distinct org_name) as org_count FROM
                                      %s
                                      WHERE org_name IN (%s)
                                      GROUP BY hfam) b
                                      ON a.hfam = b.hfam
                                      WHERE org_name IN (%s)
                                      GROUP BY a.hfam, org_name, protein_id, ipr_id;",
                                      IPRtempTable, IPR_homoTable,IPR_homoTable, organism_formattet, 
                                      organism_formattet)
        
        # Check if table already exists
        checkTable <- aspDbFetch(sprintf("SHOW TABLES LIKE '%s';", IPRtempTable))
        
        if(nrow(checkTable) < 1) {
            # If the table doesn't exist - Create table and index   
            newTable <- aspDbFetch(IPRtempTable_query)
            index <- aspDbFetch(sprintf("CREATE INDEX i_hfamProtIPR ON %s (hfam, protein_id, ipr_id);", IPRtempTable))
            index <- aspDbFetch(sprintf("CREATE INDEX i_orgProtIPR ON %s (org_name, protein_id, ipr_id);", IPRtempTable))
        } else {
            # Retrieve species names from table
            IPRtempOrg <- aspDbFetch(sprintf("SELECT DISTINCT(org_name) FROM %s;", IPRtempTable))
            IPRtempOrg <- as.vector(IPRtempOrg$org_name)
            # If some of the species are missing drop table and create new + index
            if(!(all(IPRtempOrg %in% organisms) == TRUE && all(organisms %in% IPRtempOrg) == TRUE)) {
                dropTable <- aspDbFetch(sprintf("DROP TABLE %s;", IPRtempTable))
                newTable <- aspDbFetch(IPRtempTable_query)
                index <- aspDbFetch(sprintf("CREATE INDEX i_hfamProtIPR ON %s (hfam, protein_id, ipr_id);", IPRtempTable))
                index <- aspDbFetch(sprintf("CREATE INDEX i_orgProtIPR ON %s (org_name, protein_id, ipr_id);", IPRtempTable))
            }
        }
    } else {
        sprintf('# INFO: The temporary InterPro table [%s] already exists and will be used.', IPRtempTable)
    }
    
    ###########################################################################
    # CALCULATE INTERPRO COVERAGE FOR HFAMS/PROTEINS
    ##########################################################################
    
    # -------------------------------------------------------------------------
    # For proteins: Per species - calculate ipr coverage of pan/core/unique
    # -------------------------------------------------------------------------
    # PAN PROTEINS
    pan_IPR_query <- sprintf("SELECT tc.org_name, COUNT(DISTINCT tc.org_name, tc.protein_id) as 'Pan_InterPro_annotated' 
                             FROM (                   
                             SELECT org_name, protein_id, ipr_id
                             FROM %s
                             GROUP BY org_name, protein_id
                             HAVING GROUP_CONCAT(ipr_id) IS NOT NULL) tc 
                             group by org_name;", IPRtempTable)
    
    pan_NULL_query <- sprintf("SELECT tc.org_name, COUNT(DISTINCT tc.org_name, tc.protein_id) as 'Pan_Not_annotated' 
                              FROM (                   
                              SELECT org_name, protein_id, ipr_id
                              FROM %s
                              GROUP BY org_name, protein_id
                              HAVING GROUP_CONCAT(ipr_id) IS NULL) tc 
                              group by org_name;", IPRtempTable)
    
    # UNIQUE PROTEINS
    unique_IPR_query <- sprintf("SELECT tc.org_name, COUNT(DISTINCT tc.org_name, tc.protein_id) as 'Unique_InterPro_annotated' 
                                FROM (                   
                                SELECT org_name, protein_id, ipr_id
                                FROM %s
                                WHERE org_count = 1
                                GROUP BY org_name, protein_id
                                HAVING GROUP_CONCAT(ipr_id) IS NOT NULL) tc 
                                group by org_name;", IPRtempTable)
    
    unique_NULL_query <- sprintf("SELECT tc.org_name, COUNT(DISTINCT tc.org_name, tc.protein_id) as 'Unique_Not_annotated' 
                                 FROM (                   
                                 SELECT org_name, protein_id, ipr_id
                                 FROM %s
                                 WHERE org_count = 1
                                 GROUP BY org_name, protein_id
                                 HAVING GROUP_CONCAT(ipr_id) IS NULL) tc 
                                 group by org_name;", IPRtempTable)
    
    # CORE PROTEINS
    core_IPR_query <- sprintf("SELECT tc.org_name, COUNT(DISTINCT tc.org_name, tc.protein_id) as 'Core_InterPro_annotated' 
                              FROM (                   
                              SELECT org_name, protein_id, ipr_id
                              FROM %s
                              WHERE org_count = %s
                              GROUP BY org_name, protein_id
                              HAVING GROUP_CONCAT(ipr_id) IS NOT NULL) tc 
                              group by org_name;", IPRtempTable, length(organisms))
    
    core_NULL_query <- sprintf("SELECT tc.org_name, COUNT(DISTINCT tc.org_name, tc.protein_id) as 'Core_Not_annotated' 
                               FROM (                   
                               SELECT org_name, protein_id, ipr_id
                               FROM %s
                               WHERE org_count = %s
                               GROUP BY org_name, protein_id
                               HAVING GROUP_CONCAT(ipr_id) IS NULL) tc 
                               group by org_name;", IPRtempTable, length(organisms))
    
    
    
    # -------------------------------------------------------------------------
    # For hfams: Per species - calculate ipr coverage of pan/core/unique
    # -------------------------------------------------------------------------
    # UNIQUE FAMILY
    unique_family_query <- sprintf("SELECT COUNT(DISTINCT hfam) as 'Unique'
                                   FROM %s
                                   WHERE org_count = 1;", IPRtempTable)
    # PAN FAMILY
    pan_family_query <- sprintf("SELECT COUNT(DISTINCT hfam) as 'Pan'
                                FROM %s;", IPRtempTable)
    # CORE FAMILY
    core_family_query <- sprintf("SELECT COUNT(DISTINCT hfam) as 'Core'
                                 FROM %s
                                 WHERE org_count = %s;", IPRtempTable, length(organisms))
    
    
    # -------------------------------------------------------------------------
    # Run queries for proteins and hfams
    # -------------------------------------------------------------------------
    # PAN
    pan_IPR <- aspDbFetch(pan_IPR_query)
    pan_NULL <- aspDbFetch(pan_NULL_query)
    pan_family <- aspDbFetch(pan_family_query)
    # UNIQUE
    unique_IPR <- aspDbFetch(unique_IPR_query)
    unique_NULL <- aspDbFetch(unique_NULL_query)
    unique_family <- aspDbFetch(unique_family_query)
    # CORE
    core_IPR <- aspDbFetch(core_IPR_query)
    core_NULL <- aspDbFetch(core_NULL_query)
    core_family <- aspDbFetch(core_family_query)
    
    
    
    # -------------------------------------------------------------------------
    # Fill empty values in the pan/core/unique tabels
    # -------------------------------------------------------------------------
    
    ##### AUTO FILLING EMPTY VALUES #####
    # Auto-fill empty pan/core/unique values to ensure that all the orgs are included in the plots
    tablelist <- list(pan_IPR, pan_NULL, unique_IPR, unique_NULL, core_IPR, core_NULL)
    newtablist = list()
    for(tab in tablelist){
        tab_org <- tab$org_name
        new_tab = data_frame()
        for ( org in organisms ) {
            if ( !( org %in% tab_org ) ) {
                newrow <- tab[1,] # copy existing row
                newrow$org_name <- as.character(org)
                newrow[, 2] <- as.numeric('0') 
                new_tab <- data.frame(rbind(tab, newrow))
            }
        }
        if(nrow(new_tab)==0){
            newtablist[[length(newtablist)+1]] <- tab
        } else {
            newtablist[[length(newtablist)+1]] <- new_tab
        }
    }
    
    # Extract the names from the generated auto-filled tables
    pan_IPR <- newtablist[[1]]
    pan_NULL <- newtablist[[2]]
    unique_IPR <- newtablist[[3]]
    unique_NULL <- newtablist[[4]]
    core_IPR <- newtablist[[5]]
    core_NULL <- newtablist[[6]]
    
    
    # Merge and melt (rearrange) the auto-filled tables
    # OBS: Here are two versions of the tables - one for writing to csv and one for plots
    totalCORE <- merge(core_IPR,core_NULL,by="org_name")
    totalCOREmelt <- melt(totalCORE, id=c("org_name"))
    totalPAN <- merge(pan_IPR,pan_NULL,by="org_name")
    totalPANmelt <- melt(totalPAN, id=c("org_name"))
    totalUNIQUE <- merge(unique_IPR,unique_NULL,by="org_name")
    totalUNIQUEmelt <- melt(totalUNIQUE, id=c("org_name"))
    
    ###########################################################################
    # CALCULATE OVERALL PAN/CORE/UNIQUE DISTRUBTION OF PROTEINS/HFAMS
    ##########################################################################
    
    # Calculate total unique/core/pan protein counts
    Unique <- sum(colSums(Filter(is.numeric, totalUNIQUE)))
    Pan <- sum(colSums(Filter(is.numeric, totalPAN)))
    Core <- sum(colSums(Filter(is.numeric, totalCORE)))
    
    # Collect into a table "list" and rename column
    sumupProtein <- as.data.frame(rbind(Pan,Core, Unique))
    colnames(sumupProtein) <- "Genes"
    
    # Combine total unique/core/pan families counts and rename column
    sumupFamily <- as.data.frame(cbind(pan_family, core_family, unique_family))
    sumupFamily <- t(sumupFamily)
    colnames(sumupFamily) <- "Families"
    
    # Combine the two tables
    sumup <- as.data.frame(cbind(sumupProtein, sumupFamily))
    sumup$data <- rownames(sumup)
    
    # Melt table to use in plots
    sumupMelt <- melt(sumup, id=c("data"))
    sumupMelt$data <- factor(sumupMelt$data, c("Pan", "Core", "Unique"))
    
    # Converting numbers to include dot or comma 1000-separation
    if (valueSeparation == "d"){
        sumupMelt$dot <- as.factor( format(sumupMelt$value, decimal.mark=",", big.mark=".", trim=TRUE)) #DOT
    } else {
        sumupMelt$dot <- as.factor( format(sumupMelt$value, decimal.mark=".", big.mark=",", trim=TRUE)) #COMMA
    }
    
    
    ###########################################################################
    # CALCULATE IPR COVERAGE PER SPECIES OF PAN/CORE/UNIQUE + REARRANGE TABLES
    ###########################################################################
    # CORE
    # Set order for column
    names(totalCOREmelt)[names(totalCOREmelt)=="org_name"] <- "Species"
    totalCOREmelt$Species <- factor(totalCOREmelt$Species, levels=name_ordered_tips)
    # Reorder the table
    totalCOREmelt <- totalCOREmelt %>% arrange(Species)
    # Rename table
    totalCOREmelt$Species <- with(org_table, real_name[match(totalCOREmelt$Species, name)])
    totalCOREmelt$Species <- factor(totalCOREmelt$Species, levels=unique(totalCOREmelt$Species))
    
    # PAN
    # Set order for column
    names(totalPANmelt)[names(totalPANmelt)=="org_name"] <- "Species"
    totalPANmelt$Species <- factor(totalPANmelt$Species, levels=name_ordered_tips)
    # Reorder the table
    totalPANmelt <- totalPANmelt %>% arrange(Species)
    # Rename table
    totalPANmelt$Species <- with(org_table, real_name[match(totalPANmelt$Species, name)])
    totalPANmelt$Species <- factor(totalPANmelt$Species, levels=unique(totalPANmelt$Species))
    
    
    # UNIQUE
    # Set order for column
    names(totalUNIQUEmelt)[names(totalUNIQUEmelt)=="org_name"] <- "Species"
    totalUNIQUEmelt$Species <- factor(totalUNIQUEmelt$Species, levels=name_ordered_tips)
    # Reorder the table
    totalUNIQUEmelt <- totalUNIQUEmelt %>% arrange(Species)
    # Rename table
    totalUNIQUEmelt$Species <- with(org_table, real_name[match(totalUNIQUEmelt$Species, name)])
    totalUNIQUEmelt$Species <- factor(totalUNIQUEmelt$Species, levels=unique(totalUNIQUEmelt$Species))
    
    # factorize and add dots or comma as 1000-separation to value
    if (valueSeparation == "d"){
        totalPANmelt$dot <- as.factor( format(totalPANmelt$value, decimal.mark=",", big.mark=".", trim=TRUE) ) #DOT
        totalCOREmelt$dot <- as.factor( format(totalCOREmelt$value, decimal.mark=",", big.mark=".", trim=TRUE) ) #DOT
        totalUNIQUEmelt$dot <- as.factor( format(totalUNIQUEmelt$value, decimal.mark=",", big.mark=".", trim=TRUE) ) #DOT
    } else {
        totalPANmelt$dot <- as.factor( format(totalPANmelt$value, decimal.mark=".", big.mark=",", trim=TRUE) ) #COMMA
        totalCOREmelt$dot <- as.factor( format(totalCOREmelt$value, decimal.mark=".", big.mark=",", trim=TRUE) ) #COMMA
        totalUNIQUEmelt$dot <- as.factor( format(totalUNIQUEmelt$value, decimal.mark=".", big.mark=",", trim=TRUE) ) #COMMA
    }
    
    
    
    ###########################################################################
    # PLOTTING THE FIGURES
    ###########################################################################
    print("# INFO: Creating plots")
    # -------------------------------------------------------------------------
    # Phylogram tree + Pan, core and unique bar plots incl. interPro coverage
    # -------------------------------------------------------------------------
    totalPANmelt$variable <- reorder.factor(totalPANmelt$variable, new.order = rev(unique(totalPANmelt$variable)))
    totalCOREmelt$variable <- reorder.factor(totalCOREmelt$variable, new.order = rev(unique(totalCOREmelt$variable)))
    totalUNIQUEmelt$variable <- reorder.factor(totalUNIQUEmelt$variable, new.order = rev(unique(totalUNIQUEmelt$variable)))
    
    # Phylo tree
    max_rootDist <- max(distRoot(fullTree_newNames, method = "patristic"))
    tree_width <- max_rootDist + max_rootDist/2 # calculated to 0.80 manually set to 0.76
    
    if (tree_width > 1){
        offset = 0.02
    } else {
        offset = 0.0001
    }
    
    tree <-  ggtree(fullTree_newNames) + geom_tippoint(color="black", shape=16, size=2) + geom_tiplab(offset = offset) + theme_tree2() # phylotree
    IPR_PCU_phylo <- IPRcov_PanCoreUnique_perSpecies(totalPANmelt, totalCOREmelt, totalUNIQUEmelt, tree, tree_width, plotCol)
    
    # Cladogram tree
    fullTree_brlen2 <- compute.brlen(fullTree_newNames,0.5)
    max_rootDist <- max(distRoot(fullTree_brlen2, method = "patristic"))
    tree_width <- max_rootDist + (max_rootDist/2) + 0.5 # calculated to 11.25 manually set to 12
    tree <-  ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + geom_tiplab(offset = 0.1, align = TRUE) + theme_tree2() # computed_brlen + align.tips
    IPR_PCU_cladogram <- IPRcov_PanCoreUnique_perSpecies(totalPANmelt, totalCOREmelt, totalUNIQUEmelt, tree, tree_width, plotCol)
    
    # Plot
    ggsave(paste("IPR_PanCoreUnique_perSpecies_phylo_",fileName,'.pdf', sep = ''), plot = IPR_PCU_phylo, width = 22, height = length(organisms)*0.45)
    ggsave(paste("IPR_PanCoreUnique_perSpecies_phylo_",fileName,'.svg', sep = ''), plot = IPR_PCU_phylo, width = 22, height = length(organisms)*0.45)
    
    ggsave(paste("IPR_PanCoreUnique_perSpecies_cladogram_",fileName,'.pdf', sep = ''), plot = IPR_PCU_cladogram, width = 22, height = length(organisms)*0.45)
    ggsave(paste("IPR_PanCoreUnique_perSpecies_cladogram_",fileName,'.svg', sep = ''), plot = IPR_PCU_cladogram, width = 22, height = length(organisms)*0.45)
    
    # -------------------------------------------------------------------------
    # Phylogram tree + Pan, core and unique bar plots
    # -------------------------------------------------------------------------
    # Calculate sum
    Ptotals <- totalPANmelt %>%
        group_by(Species) %>%
        summarize(total = sum(value))
    
    Ctotals <- totalCOREmelt %>%
        group_by(Species) %>%
        summarize(total = sum(value))
    
    Utotals <- totalUNIQUEmelt %>%
        group_by(Species) %>%
        summarize(total = sum(value))
    
    colnames(Ptotals) <- c("Species", "Pan")
    colnames(Ctotals) <- c("Species", "Core")
    colnames(Utotals) <- c("Species", "Unique")
    
    # Reformat sum table
    PCU_combo <- merge(as.data.frame(Ctotals), as.data.frame(Utotals), by = "Species" )
    PCU_combo <- merge(PCU_combo, as.data.frame(Ptotals), by = "Species" )
    Ptotals_all <-  PCU_combo[c("Species", "Pan")]
    PCU_combo$Accessory <- (PCU_combo$Pan - PCU_combo$Core - PCU_combo$Unique)
    PCU_combo <- PCU_combo[c("Species", "Core", "Accessory", "Unique")]
    PCU_melt <- melt(PCU_combo, id = "Species")
    PCU_melt$variable <- factor(PCU_melt$variable, levels=c( "Accessory","Unique",  "Core"))
    
    # Plot
    max_rootDist <- max(distRoot(fullTree_brlen2, method = "patristic"))
    tree_width <- max_rootDist + (max_rootDist/2) # calculated to 11.25 manually set to 12
    
    phy <-  ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + geom_tiplab(offset = 0.1, align = TRUE) + theme_tree2() # computed_brlen
    PCU_plot <- PanCoreUnique_perSpecies(PCU_melt, Ptotals, phy, plotCol, tree_width)
    
    ggsave(paste("PanCoreUnique_perSpecies_",fileName,'.pdf', sep = ''), plot = PCU_plot, width = 15, height = length(organisms)*0.45)
    ggsave(paste("PanCoreUnique_perSpecies_",fileName,'.svg', sep = ''), plot = PCU_plot, width = 15, height = length(organisms)*0.45)
    
    # -------------------------------------------------------------------------
    # Pan, core and unique bar plots total
    # -------------------------------------------------------------------------
    
    bar_dodged <- prot_hfam_PanCoreUnique(sumupMelt, xcol=sumupMelt$variable, ycol=sumupMelt$value, varcol=sumupMelt$data, x_axis_label = "Overall distribution", valueSeparation = valueSeparation, colorChoice = plotCol[c(3,1,5)])
    
    ggsave(paste("PanCoreUnique_total_",fileName,'.pdf', sep = ''), plot = bar_dodged, width = 6, height = 10)
    ggsave(paste("PanCoreUnique_total_",fileName,'.svg', sep = ''), plot = bar_dodged, width = 6, height = 10)
    
    
    
    ###########################################################################
    # TREE CALCULATIONS
    ###########################################################################
    
    ###########################################################################
    # SETTING COUNTS ON PHYLO TREE NODES
    ###########################################################################
    # Creating a list with each element being a vector of organism names
    # for each node of the full phylogeny tree
    IPRperNode_df <- data.frame(node= integer(0), IPR_fract= numeric(0), 
                                NoIPR_fract = numeric(0), IPR= integer(0), 
                                NoIPR=integer(0), totalProt = integer(0), 
                                totalClust = integer(0), orgs = character(0))
    
    fullTree_brlen2 <- compute.brlen(fullTree_newNames,0.5)
    max_rootDist <- max(distRoot(fullTree_brlen2, method = "patristic"))
    tree_width <- max_rootDist + (max_rootDist/2) + 1.5 # calculated to 11.25 manually set to 12
    p <-  ggtree(fullTree_brlen2) + geom_tiplab(offset = 1.25, align = TRUE) + theme_tree2() # computed_brlen + align.tips
    nodes <- unique(p$data$parent)
    
    for (nod in nodes){
        nodeOrgs <- extract.clade(fullTree,nod)$tip.label
        org_form <- toString(shQuote(mixedsort(nodeOrgs)))
        query_totalProtClust <- sprintf("SELECT COUNT(DISTINCT ta.hfam, org_name, protein_id) 
                                        as TotalProteins, COUNT(DISTINCT ta.hfam) as TotalClusters
                                        FROM %s AS ta
                                        JOIN (
                                        SELECT hfam
                                        FROM %s
                                        WHERE org_name IN (%s) AND org_count = %s
                                        GROUP BY hfam
                                        HAVING count(DISTINCT org_name) = %s) tb
                                        ON (ta.hfam = tb.hfam);", IPRtempTable,  IPRtempTable, org_form, length(nodeOrgs),  length(nodeOrgs))
        
        query_NULLprot <- sprintf("SELECT COUNT(DISTINCT tc.hfam, tc.org_name, tc.protein_id) as NULLproteins
                                  FROM (SELECT ta.hfam, org_name, protein_id, ipr_id
                                  FROM %s AS ta
                                  JOIN (
                                  SELECT hfam
                                  FROM %s
                                  WHERE org_name IN (%s) AND org_count = %s
                                  GROUP BY hfam
                                  HAVING count(DISTINCT org_name) = %s) tb
                                  ON (ta.hfam = tb.hfam)
                                  GROUP BY ta.hfam, protein_id
                                  HAVING GROUP_CONCAT(`ipr_id`) IS NULL) tc;", IPRtempTable,  IPRtempTable, org_form, length(nodeOrgs),  length(nodeOrgs))
        totalProtCluster <- aspDbFetch(query_totalProtClust)
        totalClust <- totalProtCluster[[2]]
        totalProt <- as.numeric(totalProtCluster[[1]])
        NULLcovProt <- as.numeric(aspDbFetch(query_NULLprot))
        IPRcovProt <- totalProt-NULLcovProt
        df <-  data.frame(node= nod, IPR_fract= IPRcovProt/totalProt, NoIPR_fract =  NULLcovProt/totalProt, IPR= IPRcovProt, NoIPR =  NULLcovProt, totalProt = totalProt, totalClust = totalClust, orgs = org_form)
        IPRperNode_df <- rbind(IPRperNode_df, df)
    }
    
    
    tipNodes <- p$data[ which(p$data$isTip=='TRUE'), ]$node
    for (tipNode in tipNodes){
        tipOrgs <- p$data[ which(p$data$node==tipNode), ]$label
        tipOrgs <- org_table[match(tipOrgs, org_table$real_name),1]
        
        query_totalProtClust <- sprintf("SELECT COUNT(DISTINCT hfam, org_name, protein_id) 
                                        as TotalProteins, COUNT(DISTINCT hfam) as TotalClusters
                                        FROM %s 
                                        WHERE org_name = '%s' AND org_count = 1;",IPRtempTable, tipOrgs)  
        
        query_NULLprot <- sprintf("SELECT COUNT(DISTINCT hfam, protein_id) as NULLproteins
                                  FROM (                                   
                                  SELECT hfam, protein_id, ipr_id
                                  FROM %s 
                                  WHERE org_name IN ('%s') AND org_count = 1
                                  GROUP BY hfam, protein_id
                                  HAVING GROUP_CONCAT(`ipr_id`) IS NULL) ta;",IPRtempTable, tipOrgs)
        
        
        totalProtCluster <- aspDbFetch(query_totalProtClust)
        totalClust <- totalProtCluster[[2]]
        totalProt <- as.numeric(totalProtCluster[[1]])
        NULLcovProt <- as.numeric(aspDbFetch(query_NULLprot))
        IPRcovProt <- totalProt-NULLcovProt
        df <-  data.frame(node= tipNode, IPR_fract= IPRcovProt/totalProt, NoIPR_fract =  NULLcovProt/totalProt, IPR= IPRcovProt, NoIPR =  NULLcovProt, totalProt = totalProt, totalClust = totalClust, orgs = tipOrgs)
        IPRperNode_df <- rbind(IPRperNode_df, df)
    }
    

    IPRperNode_df[is.na(IPRperNode_df)] <- 0
    p$data <- merge(p$data, IPRperNode_df[,1:7], by.x = "node",  by.y = "node",all = TRUE)
    
    #pFam <- p + geom_label2(aes(label=totalClust),  fill="#3e3e3fff", size = 3, color = "white") #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
    pFam1 <- p + geom_label2(aes(label=totalClust, subset = !isTip),  fill=plotCol[3], size = 3, color = "white") #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
    pFam2 <- pFam1 + geom_label2(aes(label=totalClust, subset = isTip),  fill="white", size = 3, color = "black", hjust=-.1) #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
    
    pProt1 <- p + geom_label2(aes(label=totalProt, subset = !isTip),  fill=plotCol[3], size = 3, color = "white") #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
    pProt2 <- pProt1 + geom_label2(aes(label=totalProt, subset = isTip),  fill="white", size = 3, color = "black", hjust=-.1) #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
    
    hfam_IPR_PCU_cladogram <- IPRcov_PanCoreUnique_perSpecies(totalPANmelt, totalCOREmelt, totalUNIQUEmelt, pFam2, tree_width, plotCol)
    protein_IPR_PCU_cladogram <- IPRcov_PanCoreUnique_perSpecies(totalPANmelt, totalCOREmelt, totalUNIQUEmelt, pProt2, tree_width, plotCol)
    
    ggsave(paste("hfamCount_IPR_PanCoreUnique_perSpecies_cladogram_",fileName,'.pdf', sep = ''), plot = hfam_IPR_PCU_cladogram, width = 22, height = length(organisms)*0.45)
    ggsave(paste("hfamCount_IPR_PanCoreUnique_perSpecies_cladogram_",fileName,'.svg', sep = ''), plot = hfam_IPR_PCU_cladogram, width = 22, height = length(organisms)*0.45)
    
    ggsave(paste("proteinCount_IPR_PanCoreUnique_perSpecies_cladogram_",fileName,'.pdf', sep = ''), plot = protein_IPR_PCU_cladogram, width = 22, height = length(organisms)*0.45)
    ggsave(paste("proteinCount_IPR_PanCoreUnique_perSpecies_cladogram_",fileName,'.svg', sep = ''), plot = protein_IPR_PCU_cladogram, width = 22, height = length(organisms)*0.45)
    
    # Pie and bar plots on cladogram
    #p1 <- ggtree(fullTree_brlen2) + geom_tiplab(offset = 0.2) + xlim(0, 13)
    #pies <- nodepie(IPRperNode_df[,1:7], cols=2:3, color=c(IPR_fract='#249dc5ff', NoIPR_fract='#3e3e3fff'), alpha=0.8)
    #pie_plot <- inset(p1, pies, width=.13)
    
    #bars <- nodebar(IPRperNode_df[,1:7], cols=4:5, color=c(IPR='blue', NoIPR='#3e3e3fff'))
    #bar_plot <- inset(p1, bars, width=.03, height=.04)
    
    #ggsave(paste("pie_tree_",fileName,'.pdf', sep = ''), plot = pie_plot, width = 6, height = length(organisms)*0.26)
    #ggsave(paste("pie_tree_",fileName,'.svg', sep = ''), plot = pie_plot, width = 6, height = length(organisms)*0.26)
    
    
    p2 <- ggtree(fullTree_brlen2) + geom_tiplab(offset = 1) + xlim(0, 14)
    p2$data <- merge(p2$data, IPRperNode_df[,1:7], by.x = "node",  by.y = "node",all = TRUE)
    
    pFam1_alone <- p2 + geom_label2(aes(label=totalClust, subset = !isTip),  fill=plotCol[3], size = 3, color = "white") #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
    pFam2_alone <- pFam1_alone + geom_label2(aes(label=totalClust, subset = isTip),  fill="white", size = 3, color = "black", hjust=-.1) #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
    
    pProt1_alone <- p2 + geom_label2(aes(label=totalProt, subset = !isTip),  fill=plotCol[3], size = 3, color = "white") #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
    pProt2_alone <- pProt1_alone + geom_label2(aes(label=totalProt, subset = isTip),  fill="white", size = 3, color = "black", hjust=-.1) #, alpha = 0.5)  # aes(label=totalClust, subset = !isTip),          
 
    ggsave(paste("hfamCount_cladogram_",fileName,'.pdf', sep = ''), plot = pFam2_alone, width = 7, height = length(organisms)*0.52)
    ggsave(paste("hfamCount_cladogram_",fileName,'.svg', sep = ''), plot = pFam2_alone, width = 7, height = length(organisms)*0.52)
    
    ggsave(paste("proteinCount_cladogram_",fileName,'.pdf', sep = ''), plot = pProt2_alone, width = 7, height = length(organisms)*0.52)
    ggsave(paste("proteinCount_cladogram_",fileName,'.svg', sep = ''), plot = pProt2_alone, width = 7, height = length(organisms)*0.52)
    
    #length(organisms)
    #Follow_phylo <- sum(p2$data$totalProt)
    #pan_proteome <- sum(Ptotals_all$Pan)
    #NOTFollow_phylo <- pan_proteome - Follow_phylo
    #df <-  data.frame(node= nod, IPR_fract= IPRcovProt/totalProt, NoIPR_fract =  NULLcovProt/totalProt, IPR= IPRcovProt, NoIPR =  NULLcovProt, totalProt = totalProt, totalClust = totalClust, orgs = org_form)
    #Data_followPhylo <- data.frame(Pan_proteome= integer(0), Following_phylogeny= integer(0), 
    #                            Not_following_phylogeny = integer(0), pct_Following_phylogeny= numeric(0), 
    #                            pct_Not_following_phylogeny=numeric(0))
    #Data_followPhylo <- data.frame(Pan_proteome = pan_proteome, Following_phylogeny = Follow_phylo, Not_following_phylogeny = NOTFollow_phylo, pct_Following_phylogeny = (Follow_phylo/Pan_proteome)*100, pct_Not_following_phylogeny = (Not_following_phylogeny/Pan_proteome)*100)
    
    #length(organisms)
    #Follow_phylo_hfam <- sum(p2$data$totalClust)
    #pan_proteome_hfam <- sumup$`Protein families`[1]
    #NOTFollow_phylo_hfam <- pan_proteome_hfam - Follow_phylo_hfam
    #(NOTFollow_phylo_hfam/pan_proteome_hfam)*100
    

    ###########################################################################
    # Deleting temperary R build table
    ###########################################################################
    if (identical(createIPRtables, FALSE)){
        checkTable <- aspDbFetch(sprintf("SHOW TABLES LIKE '%s';", IPRtempTable))
        if(nrow(checkTable) < 1) {
            dropTable <- aspDbFetch(sprintf("DROP TABLE %s",IPRtempTable))
        }
    }
    
    ###########################################################################
    # WRITE PAN/CORE/UNIQUE CSV TABLES
    ###########################################################################
    
    if(isTRUE(write_tables)) {
        
        # ---------------------------------------------------------------------
        # Overall pan/core/unique distribution for proteins and hfams
        # ---------------------------------------------------------------------
        #### WRITE TABLE TO CSV ####
        if (valueSeparation == "d"){
            write.table(sumup[c("data", "Genes", "Families")], file = sprintf("PanCoreUnique_total_%s.csv", fileName), sep = ";",
                        eol = "\n", na = "NA", dec = ",", row.names = FALSE,
                        col.names = TRUE, fileEncoding = "utf-8")
        } else {
            write.table(sumup[c("data", "Genes", "Families")], file = sprintf("PanCoreUnique_total_%s.csv", fileName), sep = ";",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, fileEncoding = "utf-8")
        }
        
        
        
        write.table(sumup[c("data", "Genes", "Families")], file = sprintf("PanCoreUnique_total_%s.csv", fileName), sep = ";",
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                    col.names = TRUE, fileEncoding = "utf-8")
        
        # Melt table to use in plots
        sumupMelt <- melt(sumup, id=c("data"))
        
        
        # ---------------------------------------------------------------------
        # IPR coverage of pan/core/unique for individual species
        # ---------------------------------------------------------------------
        # rename column to use it in adjToPhylogeny
        names(totalCORE)[names(totalCORE)=="org_name"] <- "Species"
        names(totalPAN)[names(totalPAN)=="org_name"] <- "Species"
        names(totalUNIQUE)[names(totalUNIQUE)=="org_name"] <- "Species"
        names(totalPAN)[names(totalPAN)=="Pan_InterPro_annotated"] <- "Genome_InterPro_annotated"
        names(totalPAN)[names(totalPAN)=="Pan_Not_annotated"] <- "Genome_Not_annotated"
        
        # merge tabels into one and rename column
        collected <- merge(totalPAN,totalCORE,by="Species", sort = FALSE)
        collected <- merge(collected,totalUNIQUE,by="Species", sort = FALSE)
        
        # Switch JGI org_name to real_species name column
        collected$Species <- with(org_table, real_name[match(collected$Species, name)])
        
        #### WRITE TABLE TO CSV ####
        if (valueSeparation == "d"){
            write.table(collected, file = sprintf("IPR_PanCoreUnique_%s.csv", fileName), sep = ";",
                        eol = "\n", na = "NA", dec = ",", row.names = FALSE,
                        col.names = TRUE, fileEncoding = "utf-8")
        } else {
            write.table(collected, file = sprintf("IPR_PanCoreUnique_%s.csv", fileName), sep = ";",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, fileEncoding = "utf-8")
        }
        
        # ---------------------------------------------------------------------
        # Pan/core/unique/accessory for individual species
        # ---------------------------------------------------------------------
        
        PCU_combo_writeOUT <- merge(PCU_combo, as.data.frame(Ptotals_all), by = "Species" )
        names(PCU_combo_writeOUT)[names(PCU_combo_writeOUT)=="Pan"] <- "Genome"
        
        if (valueSeparation == "d"){
            write.table(PCU_combo_writeOUT, file = sprintf("PanCoreUnique_%s.csv", fileName), sep = ";",
                        eol = "\n", na = "NA", dec = ",", row.names = FALSE,
                        col.names = TRUE, fileEncoding = "utf-8")
        } else {
            write.table(PCU_combo_writeOUT, file = sprintf("PanCoreUnique_%s.csv", fileName), sep = ";",
                        eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                        col.names = TRUE, fileEncoding = "utf-8")
        }
        
    }
    setwd(oldwd)
    print("The program is finished.")
}

# -------------------------------------------------------------------------
# PLOT TREES
# -------------------------------------------------------------------------
# Saving the tree files
# OBS: I can't save the trees in svg format
# CLADOGRAMS
#svg( file = paste("cladogram_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames, branch.length="none") + geom_tippoint(color="black", shape=16, size=2) +geom_tiplab(offset = 0.5) + xlim(0, 35) 
#dev.off()
#
#pdf( file = paste("cladogram_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames, branch.length="none") + geom_tippoint(color="black", shape=16, size=2) +geom_tiplab(offset = 0.5) + xlim(0, 35) 
#dev.off()
#
#svg( file = paste("cladogram_NOdot_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames, branch.length="none") +geom_tiplab(offset = 0.15) + xlim(0, 35) 
#dev.off()
#
#pdf( file = paste("cladogram_NOdot_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames, branch.length="none")  +geom_tiplab(offset = 0.15) + xlim(0, 35) 
#dev.off()
#
#svg( file = paste("cladogram_bootstrap_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames, branch.length="none") + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(offset = 0.5) + xlim(0, 35) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) > 99), 
#                                                         size = 3, color = "red", vjust = -0.2, hjust = 1.1) + geom_text2(aes(label=label, 
#                                                                                                                              subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) < 99), size = 3, color = "black", vjust = -0.2, hjust = 1.1)
#dev.off()
#
#pdf( file = paste("cladogram_bootstrap_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames, branch.length="none") + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(offset = 0.5) + xlim(0, 35) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) > 99), 
#                                                         size = 3, color = "red", vjust = -0.2, hjust = 1.1) + geom_text2(aes(label=label,                                                                                                                                 subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) < 99), size = 3, color = "black", vjust = -0.2, hjust = 1.1)
#dev.off()
#
## PHYLOTREES
#svg(file = paste("phylotree_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames) + geom_tippoint(color="black", shape=16, size=2) + geom_tiplab(offset = 0.02) + xlim(0,3)
#dev.off()
#
#pdf(file = paste("phylotree_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames) + geom_tippoint(color="black", shape=16, size=2) + geom_tiplab(offset = 0.02) + xlim(0,3)
#dev.off()  
#
#svg( file = paste("phylotree_NOdot_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames) + geom_tiplab(offset = 0.01) + xlim(0,3) + theme_tree2()
#dev.off()
#
#pdf( file = paste("phylotree_NOdot_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames) + geom_tiplab(offset = 0.01) + xlim(0,3) + theme_tree2()
#dev.off()  
#
#svg( file = paste("phylotree_alignTips_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames) + geom_tippoint(color="black", shape=16, size=2) + geom_tiplab(align=TRUE, linesize = 0.5, offset = 0.05) + xlim(0,3) 
#dev.off()
#
#pdf( file = paste("phylotree_alignTips_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_newNames) + geom_tippoint(color="black", shape=16, size=2) + geom_tiplab(align=TRUE, linesize = 0.5, offset = 0.05) + xlim(0,3) 
#dev.off()  
#
## CLADOGRAMS COMPUTED BRANCH LENGTH
#fullTree_brlen2 <- compute.brlen(fullTree_newNames,0.5)
#
#svg( file = paste("cladogram_brlen_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(offset = 0.1) + xlim(0, 14)
#dev.off()
#
#pdf( file = paste("cladogram_brlen_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(offset = 0.1) + xlim(0, 14)
#dev.off()
#
#svg( file = paste("cladogram_brlen_alignTips_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(align=TRUE, linesize = 0.5, offset = 0.05) + xlim(0, 14)
#dev.off()
#
#pdf( file = paste("cladogram_brlen_alignTips_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(align=TRUE, linesize = 0.5, offset = 0.05) + xlim(0, 14) 
#dev.off()
#
#svg( file = paste("cladogram_brlen_bootstrap_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(offset = 0.1) + xlim(0, 14) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) > 99), 
#                                                         size = 3, color = "red", vjust = -0.2, hjust = 1.1) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & 
#                                                                                                                                  as.numeric(label) < 99), size = 3, color = "black", vjust = -0.2, hjust = 1.1)
#dev.off()
#
#pdf( file = paste("cladogram_brlen_bootstrap_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(offset = 0.1) + xlim(0, 14) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) > 99), 
#                                                         size = 3, color = "red", vjust = -0.2, hjust = 1.1) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & 
#                                                                                                                                  as.numeric(label) < 99), size = 3, color = "black", vjust = -0.2, hjust = 1.1)
#dev.off()
#
#svg( file = paste("cladogram_brlen_alignTips_bootstrap_",fileName,'.svg', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(align=TRUE, linesize = 0.5, offset = 0.05) + xlim(0, 14) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) > 99), 
#                                                                                      size = 3, color = "red", vjust = -0.2, hjust = 1.1) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & 
#                                                                                                                                                               as.numeric(label) < 99), size = 3, color = "black", vjust = -0.2, hjust = 1.1)
#dev.off()
#
#pdf( file = paste("cladogram_brlen_alignTips_bootstrap_",fileName,'.pdf', sep = ''), width = 7, height = 20) 
#ggtree(fullTree_brlen2) + geom_tippoint(color="black", shape=16, size=2) + 
#    geom_tiplab(align=TRUE, linesize = 0.5, offset = 0.05) + xlim(0, 14) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) > 99), 
#                                                                                      size = 3, color = "red", vjust = -0.2, hjust = 1.1) + geom_text2(aes(label=label, subset = !isTip & !is.na(as.numeric(label)) & 
#                                                                                                                                                               as.numeric(label) < 99), size = 3, color = "black", vjust = -0.2, hjust = 1.1)
#dev.off()
