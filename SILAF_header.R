# External script encoding functions used throughout the SILAF code.

# Mean functions cannot cope with NA values in apply, therefore redefine the mean function.
## @knitr mean.na
mean.na <- function(x){
  mean(x, na.rm = T)
}

# Normalizing data by hand-function.
## @knitr normalize
normalize <- function(x){
  sums <- apply(x, 2, function(y){sum(y, na.rm = TRUE)})
  factor <- mean(sums)/sums
  x <- data.frame(t(t(x)*factor))
}

# Define a common mode to choose experimental data >>>
## @knitr read.ms
read.ms <- function(data = total, meta = info, name = "", type = "SILAF", inf.clearing = "FALSE", sum.columns = "TRUE", 
                    normalize = "FALSE", plots = "FALSE", stats, stat.test = "X.Log.t.test.p.value", diff.value = "t.test.difference"){
  
  # Usage of the function, for instance:
  # read.ms(data = total, meta = info, name = "QC.SILAF.MONSTER", type = "SILAF", inf.clearing = "to.NA", sum.column = "FALSE",
  # plots = "TRUE", stats = monster.SILAF.stats, normalize = "FALSE")
  
  # Specifying "sub" depending on whether it is an LFQ or SILAF experiment
  if(type == "SILAF"){
    sub <- "Intensity.*\\."
  }
  if(type == "LFQ"){
    sub <- "LFQ.intensity.*\\."
  }
  
  # Which files correspond to this experiment? >
  files.current <- meta[meta$Experiment==name,][,"Tube"]
  colnames <- sub(sub, "", colnames(data)) #<Obtain columns that contain the quantification data.>
  dataset.current <- data[,colnames %in% files.current]
  
  # Eliminate all rows with only null values >
  dataset.current <- dataset.current[apply(dataset.current, 1, sum) != 0,] 
  dataset.initial <- dataset.current
  
  ## Deletes the MaxQuant calculated sum columns of light + heavy. These can be pretty annoying later. >>
  if(sum.columns == "FALSE"){
    exclude <- sub(".*\\..\\..*", "", colnames(dataset.current)) != ""
    dataset.current <- dataset.current[,exclude=="FALSE"]
    if(normalize == "TRUE"){
      dataset.current <- normalize(dataset.current)
    }
  }
  
  ## Specifying how to deal with inf values caused by log2 transformations >>
  # Leave the inf values as they are >
  if(inf.clearing == "FALSE"){ 
    dataset.out <- log2(dataset.current) # Log2 transformation
  }
  
  # Enables mean intensity calculation without including values below the detection limit >
  if(inf.clearing == "to.NA"){ 
    dataset.out <- do.call(data.frame, lapply(log2(dataset.current), 
                                              function(x){replace(x,is.infinite(x), NA)})) 
    # < Replace infinite values by NA and convert back to a data.frame
    rownames(dataset.out) <- rownames(dataset.current)
  }
  
  # Convertes inf values to 0. Applicable if we want to keep 0s as 0s >
  if(inf.clearing == "to.0"){ 
    dataset.out <- do.call(data.frame, lapply(log2(dataset.current), 
                                              function(x){replace(x,is.infinite(x), 0)})) 
    # < Replace infinite values by 0 and convert back to a data.frame
    rownames(dataset.out) <- rownames(dataset.current)
  }
  
  ## Stands independend of inf-fates due to the "else" condition at the end and prepares a list for an MA plot 
  ## option >>
  if(plots == "TRUE"){ 
    
    # Step 1: exclusion of sum light+heavy columns >
    exclude <- sub(".*\\..\\..*", "", colnames(dataset.initial)) != ""
    dataset.current <- dataset.initial[,exclude=="FALSE"]
    
    # Step 2/1: replace inf values caused by log2 transformation to NA, required for sum.na >
    dataset.current.df <- do.call(data.frame, lapply(log2(dataset.current), 
                                                     function(x){replace(x,is.infinite(x), NA)})) 
    rownames(dataset.current.df) <- rownames(dataset.current)
    dataset.current <- dataset.current.df
    
    # Step 2/2: Label swap correction if the replicates indicate an inverse labelling
    # Which files have "inv" in their name?
    #inv.files <- intersect(files.current, meta[sub(".*inv","", meta[,"Replicate"]) == "",][,"Tube"])
    #Corresponding files in current dataset
    #inv.current <- sub(".*\\.", "", colnames(dataset.current)) %in% inv.files
    # Make a logical string with +/-1
    #inv.current <- replace(inv.current, inv.current == TRUE, -1)
    #inv.current <- replace(inv.current, inv.current == 0, 1)
    # Manipulate current dataset
    #dataset.current <- data.frame(t(inv.current * t(dataset.current)))
    
    # Step 3: Include only those proteins that are also present in the Perseus-generated stats file
    overlap.current <- which(rownames(dataset.current) %in% stats[,"Protein.IDs"])
    dataset.current <- dataset.current[overlap.current,]
    overlap.current <- which(stats[,"Protein.IDs"] %in% rownames(dataset.current))
    stats <- stats[overlap.current,]
    
    # Step 4: Normalize data
    dataset.current <- normalize(dataset.current)
    
    # Step 5: Calculate the mean intensity value >
    mean.current <- apply(dataset.current, 1, mean.na)
    
    # Step 6: Add logfc to plots list element
    
    logfc.current <- apply(stats[,1:length(files.current)], 1, mean.na)
    
    ## Which column contains the difference values?
    diff.value <- which(grepl("Diff", colnames(stats)), arr.ind = TRUE)
    
    plots.current <- data.frame(mean.int = mean.current, logfc = logfc.current, padj = stats[,stat.test], Gene.names = stats[,"Gene.names"],
                                t.difference = stats[,diff.value])
    rownames(plots.current) <- rownames(dataset.current)
    
    # Step 7: Make a list with [[int]] the raw intensity values and [[means]] a named means vector >
    list.current <- list()
    list.current[["int"]] <- dataset.out
    list.current[["plots"]] <- plots.current
    list.current[["stats"]] <- stats
    
    # Step 6: Return the list ONLY IF plots == "TRUE" >
    return(list.current)
  } else {
    return(dataset.out)} # < ... or a simple data.frame
}

# Scatterplot function. Requires sum.columns = "TRUE" >>>
## @knitr scatter.SILAF
scatter.SILAF <- function(data, meta){
  
  ## Prepare the data >>
  file.current <- sub(".*\\.","",colnames(data)[1]) # What is the filenumber?
  R.current <- abs(cor(data[,2], data[,3], use = "complete.obs")) # Correlation coefficient R
  
  print(paste("R of", file.current, "is:", round(R.current, digits = 3))) 
  # < Some feedback keeps the mood up
  
  ## Plotting >>
  ggplot(data = data, aes(x = log2(data[,1]), y = log2(data[,2]))) +
    geom_point() +
    geom_density2d(aes(colour=..level..)) + # Scatterplot
    #stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 500) +
    #scale_fill_continuous(low = color.transparent, high = "#073878") +
    #geom_point(alpha = 0.01, shape = 20) +
    labs(title = paste("L/H of Replicate", meta[meta[,"Tube"] == file.current,]$Replicate), 
         x = colnames(data)[1], y = colnames(data)[3]) +
    scale_colour_gradient(low="white",high="orange") +
    geom_abline(slope = 1, intercept = 0, color = "grey") +
    theme_classic() +
    theme(legend.position="none") +
    annotate("text", x = 25, y = 35, 
             label = paste("italic(R) ==", round(R.current, digits = 4)), parse = TRUE)+
    annotate("text", x = 30, y = 20, 
             label = paste("italic(N0) ==", sum(data[1] == data[3], na.rm = T)), parse = TRUE)+
    annotate("text", x = 20, y = 30, 
             label = paste("italic(N0) ==", sum(data[1] == data[2], na.rm = T)), parse = TRUE, 
             angle = 90) +
    xlim(20,37)+
    ylim(20,37)+
    coord_fixed()
}

# Volcano plot function. Requires the datalist >>>
## @knitr volcano.plot
volcano.plot <- function(data, x = "logfc", y = "padj", title = NULL, xintercept = 1, sig.column, threshold = 1){
  color <- data[["stats"]][,sig.column]=="+" & abs(data[["plots"]][,"logfc"]) > threshold   # < for color scale of significant hits
  color <- replace(color, is.na(color) == TRUE, "FALSE")
  data <- data[["plots"]]
  
  
  ggplot(data = data, aes(x = data[,x], y = data[,y])) +
    geom_point(aes(color = color)) +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(
      data = data[color == TRUE,],
      aes(x = (data[color == TRUE,])[,x], y = (data[color == TRUE,])[,y], 
          label = (data[color == TRUE,])[,"Gene.names"]),
      size = 4,
      box.padding = 0.25,
      point.padding = 0.3) +
    theme_bw() +
    labs(title = title, x = "log2(fold change)", y = "-log10(adjusted p-value)") +
    theme(legend.position = "none") +
    geom_vline(aes(xintercept = -xintercept), lty = "dotted") +
    geom_vline(aes(xintercept = xintercept), lty = "dotted")+
    coord_fixed(ratio = 1)
}

# Histogram function >>>
panel.hist.add <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE) # set breaks here >
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

# Histogram function >>>
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = TRUE) # set breaks here >
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

# Smooth Scatterplot function, add TRUE >>>
panel.smooth.add <- function(x, y){
  smoothScatter(x,y, add=TRUE)
  r <- round(cor(x, y, use = "complete.obs"), digits=4)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.4, 0.85, txt)
  ylim=c(20,40)
  xlim=c(20,40)
  abline(0,1, col = "grey")
}

# Smooth Scatterplot function, add FALSE >>>
panel.smooth <- function(x, y){
  smoothScatter(x,y, add=FALSE)
  r <- round(cor(x, y, use = "complete.obs"), digits=4)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.4, 0.85, txt)
  ylim=c(20,40)
  xlim=c(20,40)
  abline(0,1, col = "grey")
}

# MA plot, requires the datalist >>>
plot.MA <- function(data, sig.column = "t.test.Significant", title = NULL, threshold = 1){
  # For color scale of significant hits:
  color <- data[["stats"]][,sig.column]=="+" & abs(data[["plots"]][,"logfc"]) > threshold   
  color <- replace(color, is.na(color) == TRUE, "FALSE")
  
  ggplot(data = data[["plots"]], aes(x = log2(mean.int), y = logfc))+
    geom_point(aes(color = color))+
    scale_color_manual(values = c("grey", "red"))+
    labs(title = title, x= "log2 intensity", y = "log2 fold change")+
    theme_bw()+
    geom_hline(aes(yintercept = 0), lty = "dotted")+
    geom_hline(aes(yintercept = +1), lty = "dotted")+
    geom_hline(aes(yintercept = -1), lty = "dotted")+
    theme(legend.position = "none")
}

# GSEA data table
## @knitr gsea.table
gsea.table <- function(data.list, name, column = ""){
  logfc.table <- data.frame(ID = sub(";.*","",rownames(data.list[["plots"]])), 
                            logfc = data.list[["plots"]][,"logfc"])
  write.table(logfc.table, file = name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Redundant function using the writing format of gsea.table
## @knitr ipa
ipa <- function(x, name){
  write.table(x, file = paste("IPA.", name, sep = ""), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}

# LFQ log2 fold change calculation
## @knitr lfq.qc
lfq.fc <- function(data, meta, ctrl, treat, stats,plots = "FALSE"){
  # Narrow down to proteins that are in the Preseus stats file
  data.current <- data[rownames(data) %in% stats[,"Majority.protein.IDs"],] 
  
  # Returns the ones that are missing in the stats files. Should be mainly  contaminants
  print(paste("Missing:", (stats[(stats[,"Majority.protein.IDs"] %in% rownames(data.current))==FALSE,])
              [,"Majority.protein.IDs"])) 
  
  # Narrow down the stats file by the missing ones
  stats <- stats[stats[,"Majority.protein.IDs"] %in% rownames(data.current),]
  
  # Reduce the meta file
  meta.current <- meta[(meta[,"Tube"] %in% as.numeric(sub("LFQ\\.intensity\\.", "", colnames(data.current)))) 
                       & !is.na(meta[,"Condition"]),] 
  
  # Replace 0s by NAs to calculate proper mean values:
  data.current <- do.call(data.frame, lapply(data.current, function(x){replace(x,x == 0, NA)}))
  
  # Calculate mean values for each condition:
  mean.ctrl <- apply(data.current[,(meta.current[,"Condition"] %in% ctrl)], 1, mean.na)
  mean.treat <- apply(data.current[,(meta.current[,"Condition"] %in% treat)], 1, mean.na)
  mean.total <- apply(data.current, 1, mean.na)
  
  # Calculate the log2 fold change:
  logfc.current <- log2(mean.treat/mean.ctrl)
  
  if(plots == "FALSE"){return(mean.current)}
  if(plots == "TRUE"){
    plots.current <- data.frame(mean.int = mean.total, logfc = logfc.current, 
                                padj = stats[,"X.Log.welch.p.value"], Gene.names = stats[,"Gene.names"])
    rownames(plots.current) <- rownames(data.current)
    
    list.current <- list()
    list.current[["int"]] <- data.current
    list.current[["plots"]] <- plots.current
    list.current[["stats"]] <- stats
    
    return(list.current)
  }
}

fraction.identifier <- function(rep1, rep1.name, rep2, rep2.name, folder, fraction.top = "FALSE", fraction.low = "FALSE"){
  
  if(fraction.top == "FALSE"){
    fraction.top <- readline("What is upper cutoff? > ")
    fraction.low <- readline("What is lower cutoff? > ")
  }
  
  t.flies.rep1 <- rep1[(rep1$propheavy>=fraction.low & rep1$propheavy<=fraction.top),] # Narrow down to proteins OI in #rep1
  t.flies.rep2 <- rep2[(rep2$propheavy>=fraction.low & rep2$propheavy<=fraction.top),] # Narrow down to proteins OI in #rep2
  t.flies.temp <- rownames(t.flies.rep1[rownames(t.flies.rep2) %in% rownames(t.flies.rep1),]) # Narrow down to proteins in both samples below 0.9
  t.flies.temp <- sub(";.*", "", t.flies.temp) # Take the first part of the protein identifier
  t.flies.temp <- select(up, t.flies.temp, "ENTREZ_GENE")
  t.flies.temp <- t.flies.temp$ENTREZ_GENE
  t.flies.temp <- t.flies.temp[!(sub("\\..*", "", t.flies.temp) == "NA")] # Removes a couple of NA values at the end of the list
  
  
  ## Reference geneset:
  t.flies.temp.reference <- c(rownames(rep1), rownames(rep2)) # Combine the protein list of both replicates
  t.flies.temp.reference <- unique(t.flies.temp.reference) # Narrow down to the ones in common
  t.flies.temp.reference <- sub(";.*", "", t.flies.temp.reference) # Take the first part of the protein identifier
  t.flies.temp.reference <- select(up, t.flies.temp.reference, "ENTREZ_GENE")
  t.flies.temp.reference <- t.flies.temp.reference$ENTREZ_GENE
  t.flies.temp.reference <- t.flies.temp.reference[!(sub("\\..*", "", t.flies.temp.reference) == "NA")] # Removes a couple of NA values at the end of the list
  
  ## Write files:
  write.table(t.flies.temp, file = paste(folder, "/", rep1.name, ".", rep2.name, ".", fraction.low, ".", fraction.top, ".txt", sep = ""), col.names = F, row.names = F, quote = F)
  write.table(t.flies.temp.reference, file = paste(folder, "/", rep1.name, ".", rep2.name, ".", fraction.low, ".", fraction.top, "reference.txt", sep = ""), col.names = F, row.names = F, quote = F)
}

# Mapping of KEGG members in a labelling efficiency to intensity plot
kegg.mapping <- function(data, suffix, pathway, species = "dme"){
  data.entrez <- select(up, sub(";.*", "", rownames(data)), "ENTREZ_GENE")
  data$entrez <- data.entrez[match(sub(";.*", "", rownames(data)), data.entrez$UNIPROTKB),]$ENTREZ_GENE
  
  pathview.tmp <- data$propheavy
  names(pathview.tmp) <- data$entrez
  
  kegg.id <- keggConv(species, "up")
  
  pathway.out <- pathview(gene.data =  pathview.tmp, pathway.id = pathway,
                          species = species, gene.idtype="entrez", out.suffix = suffix, 
                          limit=list(gene=5), low="blue", high="red", mid="gray", bins=list(gene=50))
  
  pathway.mb.tmp <- sub("up:","", names(kegg.id[match(pathway.out$plot.data.gene$kegg.names, 
                                                      sub(paste(species, ":", sep = ""), "", kegg.id))]))
  
  pathway.mb.tmp.entrez <- select(up, pathway.mb.tmp, "ENTREZ_GENE")
  
  interest.tmp <- data[data$entrez %in% pathway.mb.tmp.entrez$ENTREZ_GENE,]
  
  plot <- ggplot(data = data, aes (x = propheavy, y = log10(data[,1]))) +
    geom_point(alpha = 0.2) +
    geom_point(data = interest.tmp, aes (x = propheavy, y = log10(interest.tmp[,1])), color = "firebrick3")+
    geom_vline(xintercept = mean.na(interest.tmp$propheavy), color = "firebrick3")+
    theme_classic()+
    labs(title = paste(suffix), x = "Labelling efficiency", "log10(Intensity)")
  
  return(plot)
}

# Make SD function compatible with NA values to be used in an apply-context
sd.na <- function(x){
  sd(x, na.rm = TRUE)
}