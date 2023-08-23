DEG_edgeR_func_comb <- function(out_name=NULL,folder_name,pval,logFC,out_dir,met_dir,plot_MDS,numCores,group_vect){
  
  ##############################################################################
  #Loading libraries
  library(ggplot2);library(edgeR);
  library(tximport);library(parallel)
  ##############################################################################
  # Detecting number of cores for parallelizing
  x_det <- NULL
  x_det <- parallel::detectCores()
  message(paste("Detecting cores, total cores are: ",x_det))
  
  ##############################################################################
  if (numCores > x_det) {
    stop("Number of cores exceed the maximum allowed by the machine,
         use a coherent number of cores such as 4")
  }
  
  ##############################################################################
  #Creating output dir folder
  out_dir_1 <- paste0(out_dir,"/",folder_name)
  if(!dir.exists(out_dir_1)){
    dir.create(paste0(out_dir,"/",folder_name))
  }

  # Defining plot dirs for contrasts
  plt_dir <- paste0(out_dir_1,"/graphics")
  if(!dir.exists(plt_dir)){
    dir.create(plt_dir)
  }
  
  # Defining csv dir for contrasts
  csv_dir <- paste0(out_dir_1,"/csv")
  if(!dir.exists(csv_dir)){
    dir.create(csv_dir)
  }
  
  ##############################################################################
  #Reading salmon files and combining salmon_tx2gene.tsv files
  message("Reading salmon files")

  tx2gene <- list()
  for(i in 1:length(folder_name)){
    tx2gene[[i]] <- read.table(paste0(data_dir[[i]],"/","salmon_tx2gene.tsv"))
  };rm(i)
  tx2gene <- do.call(cbind,tx2gene)
  colnames(tx2gene) <- paste("V",1:ncol(tx2gene))
  
  #reading metadata
  #tx2gene
  meta <- list()
  for(i in 1:length(folder_name)){
    meta[[i]] <- read.csv(paste0(met_dir,"/","multiqc_salmon.txt"), header = TRUE, sep="\t")
    #meta[[i]]$folder_name <- folder_name[[i]]
  };rm(i)
  
  meta <- do.call(rbind,meta)
  meta$trt <- sapply(strsplit(meta$Sample,"_"),"[[",1)
  
  # check if controls are used to compare vs. all drought treatments
  control_groups <- grep("CONTROL", meta$trt, value = TRUE)
  drought_groups <- grep("DROUGHT", meta$trt, value = TRUE)
  
  # variables that might change depending on if the samples have time or intensity
  has_time = FALSE
  has_intensity = FALSE
  
  if (length(control_groups) != length(drought_groups)) {
    # check if it's because of intensity of treatment or because of time-points
    drought_fields <- c(length(drought_groups))
    for (i in seq_along(drought_groups)) {
      sample <- drought_groups[i]
      sample_fields <- unlist(strsplit(paste(sample, ".", sep = ""), ".", fixed = TRUE))
      drought_fields[i] <- list(sample_fields)
    has_time <- any(sapply(drought_fields, function(i) i[[6]] != ""))
    has_intensity <- any(sapply(drought_fields, function(i) i[[3]] != ""))
  }}
  
  ########### FILTERING BY MAPPING RATES #######################################
  # outliers if the mapping rate is < 70% -> remove
  # not outliers if the mapping rate is < 65% -> remove
  mapping_rates <- meta$percent_mapped
  summary_rates <- summary(mapping_rates)
  q1 <- summary_rates["1st Qu."][[1]]
  q3 <- summary_rates["3rd Qu."][[1]]
  iqr <- q3-q1
  min_value <- q1-1.5*iqr # get the minimum value of the boxplot
  meta <- meta[meta$percent_mapped > min_value & meta$percent_mapped > 65 | meta$percent_mapped > 70,]
  ##############################################################################
  
  ##### AFTER FILTERING REMOVE SAMPLES THAT DO NOT HAVE DROUGHT VS CONTROL #####
  trt_unique <- unique(meta$trt)
  
  control_groups <- grep("CONTROL", trt_unique, value = TRUE)
  drought_groups <- grep("DROUGHT", trt_unique, value = TRUE)
  
  control_stage <- gsub("\\.CONTROL", ".", control_groups)
  drought_stage <- gsub("\\.DROUGHT", ".", drought_groups)
  
  if (has_time){
    drought_stage <- gsub("\\.[^.]+$", ".", drought_stage)
  }
  
  if (has_intensity){
    for (i in length(drought_stage)) {
      new_stage <- unlist(strsplit(paste(drought_stage[i], ".", sep = ""), "\\."))
      new_stage[3] <- ""
      drought_stage[i] <- paste(new_stage, collapse = ".")
    }
  }
  
  common_groups <- intersect(control_stage, drought_stage) # USE THIS TO FILTER (WHAT TO DO WHEN IS SEVEREDROUGHT...?)
  
  # take just the meta rows whose samples are in the common groups
  samples <- sapply(strsplit(meta$Sample,"_"),"[[",1)
  matches <- logical(length(meta$Sample))
  
  for (i in seq_along(samples)) {
    sample <- samples[i]
    sample_fields <- unlist(strsplit((paste(sample, ".", sep = "")), ".", fixed = TRUE))
    
    # if all controls are used to compare vs. all drought treatments (remove the intensity of the treatment)
    if (has_time) {
      sample_fields[6] = ""
    }
    
    if (has_intensity) {
      sample_fields[3] = ""
    }
    sample_fields[4] = ''
    sample_compare <- paste(sample_fields, collapse = ".")
    matches[i] <- any(sample_compare == common_groups)
  }
  
  meta <- meta[matches,]
  
  ##############################################################################
  
  samples <- sapply(strsplit(meta$Sample,"_"),"[[",2)

  #reading quant files
  list1 = list.files(path = data_dir,pattern ="quant.sf",recursive = T,full.names = T)
  list1 = list1[as.numeric(row.names(meta))]
  #original approach
  #list1 = paste0(data_dir,"/", meta$Sample,"/", "quant.sf")

  #read salmon tximport
  txi <- tximport::tximport(files=list1, type = "salmon", tx2gene = tx2gene)
  data <- as.data.frame(txi$counts)
  ab <-  as.data.frame(txi$abundance)
  df <- colSums(data)
  
  write.table(df, paste0(out_dir_1,"/","library_size.txt"), col.names = F, quote = F)
  
  ##############################################################################
  #Testing that groups have more than one sample. (this avoid NA dispersion values)
  trt_summary <- as.data.frame(tapply(meta$trt,meta$trt,length))
  trt_summary$group <- row.names(trt_summary)
  trt_summary <- trt_summary[,c(2,1)]
  colnames(trt_summary) <- c("group","count")
  
  #Using group file if it is available
  if(!is.null(group_vect)){
    if(length(meta$trt)!=length(group_vect)){
      stop("You are providing a group_vect object that not match with the number of samples.
           Please cheack and resubmit")
    } else {
      meta$trt <- group_vect
    }
  
  ##### IF WE FILTER BY MAPPING RATES THERE WILL BE TREATMENTS WITH ONLY 1  ####
  ############## SAMPLE, SO BY NOW I DECIDED TO DROP THIS OPTION ###############
  # } else {
  #   #If a group have only one sample the function stops!
  #   if(any(trt_summary$count==1)){
  #     stop("Each group must have at least two samples, please check and provide
  #          an object group_vect with the group and run again.")
  #   }
  }
  
  print(dim(data))
  print(dim(meta))
  
  ##############################################################################
  message("Defining groups with the metadata provided and levels")
  
  # Defining treatments
  colnames(data) <- meta$Sample
  colnames(ab) <- meta$Sample
  trt <- factor(meta$trt)#,levels = levels)
  
  #Defining possible combinations
  cmb_un <- as.data.frame(t(combn(levels(trt),2)))
  colnames(cmb_un) <- c("SOURCE","TARGET")
  cmb_un$COMP <- paste0(cmb_un$TARGET, "-", cmb_un$SOURCE)
  
  message("Possible combinations available: ",nrow(cmb_un))
  print(cmb_un)
  # Save normalized (but not filtered) CPM
  y <- edgeR::DGEList(counts = data,group = trt)
  y$samples$lib.size <- colSums(y$counts)
  y <- edgeR::calcNormFactors(y) # normally, norm.factors are all 1, but with this
                                 # they depend on the library size
  
  message("Saving raw data obtained (tpm and cpm)")
  
  #saving cpm
  write.table(cpm(y), paste0(out_dir_1,"/","CPM_normalized.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  #save tpm
  write.table(ab, paste0(out_dir_1,"/","abundance_drought_tpm_genelevel.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  
  ##############################################################################
  # Filtering
  message("Filtering using cpm 2")
  keep <- rowSums(cpm(data)>1)>=2 
  kept_data <- data[keep,]
  ##############################################################################
  #normal factors
  message("Normalizing")
  # Normalization and MDS plot
  d <- edgeR::DGEList(kept_data,group =trt)
  d$samples$lib.size <- colSums(d$counts)
  d <- edgeR::calcNormFactors(d)
  
  ##############################################################################
  #Plotting MDS
  t <- limma::plotMDS(d,plot = F)
  MDS_df <- data.frame(t$x, t$y)
  
  fields <- strsplit(meta$trt, '\\.')
  treatment <- unlist(lapply(fields, function(x) x[4]))  
  Features <- unlist(lapply(fields, function(x) paste(x[-c(4, which(x == ""))], collapse = ".")))
  
  #saving MDS with shapes (it plots numbers)
  
  MDS<- ggplot(data = MDS_df, aes(x = t.x, y = t.y, color = Features, shape = treatment, label = samples)) +
    geom_point() +
    geom_text(size = 3, color = "black", hjust = 0.9, vjust = -0.8) +
    theme_bw() +
    #    aes(label = meta$Sample)+ #allpoints) +
    theme(legend.box = "horizontal") +
    xlab("Leading logFC dim1") +
    ylab("Leading logFC dim2") +
    ggtitle("") +
    theme(panel.grid = element_blank())
  
  MDS <- MDS +
    guides(shape = guide_legend(title = "Treatment"),
           color = guide_legend(title = "Features"))
  
  ggsave(paste0(out_dir_1,"/","MDS.pdf"), MDS,height = 8, width = 8)
  
  if(plot_MDS==TRUE){
    plot(MDS)  
  }
  ##############################################################################
  #Refining design to create contrasts
  #design
  trt <- trt
  #using treatments and 0 to do easily
  design <- model.matrix( ~0 + trt )
  colnames( design ) <- levels( trt )
  
  ##############################################################################
  message("Estimate dispersion")
  #estimateDisp
  d <- edgeR::estimateDisp(y = d, design = design,robust = T) 
  # this function finds the best dispersion parameters that explain the observed count data
  pdf(paste0(out_dir_1,"/","plotBCV.pdf"),height = 8, width = 8)
  plotBCV(d)
  dev.off() 
  
  ##############################################################################
  #contrasts
  message("Creating contrasts")
  
  message("Using automatical combinations, analyze carefully!")
  if(nrow(cmb_un)>1){
    contrasts <- list()
    for(i in 1:nrow(cmb_un)){
      contrasts[[i]] <- makeContrasts(contrasts =  cmb_un$COMP[[i]],
                                      levels=levels(trt))
    };rm(i)
  } else {
    contrasts <- list(makeContrasts( contrasts = cmb_un$COMP[[1]] ,
                                     levels=levels(trt)))
  }
  
  for (i in length(trt)){
    message(paste("Contrast: "), trt[[i]])
  }
  
  ##############################################################################
  #obtaining DEG per contrast 
  message(paste0("Calculating DE genes in parallel, using ",numCores," cores"))
  
  if(length(contrasts)==1){
    warning("Only one contrast will be run")
  }
  
  if(length(contrasts)>0){
    
    #saving contrasts design
    x_cont <- do.call(cbind,contrasts)
    write.csv(x_cont, file= paste0(out_dir_1,"/","contrasts.csv"),na = "",row.names = T)
    
    
    #Using parallel approach
    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist=c("contrasts","d","design",
                                          "pval","logFC","cmb_un","plt_dir","out_dir_1",
                                          "csv_dir"),envir=environment())
    
    cont_results <- parallel::parLapplyLB(cl,
                                          X = seq_len(length(contrasts)),
                                          fun = function (i){
                                            #glmQLFit
                                            fit <- edgeR::glmQLFit(y = d, design = design,contrast=contrasts[[i]]) 
                                            pdf(paste0(plt_dir,"/",cmb_un$COMP[[i]],"_BCV.pdf"),height = 8, width = 8)
                                            edgeR::plotQLDisp(fit)
                                            dev.off()
                                            
                                            #testing glmQFTest
                                            qlf.2vs1 <- edgeR::glmQLFTest(glmfit = fit,contrast = contrasts[[i]])
                                            qlf.2vs1$table$FDR <- p.adjust(qlf.2vs1$table$PValue,method = "BH")
                                            qlf.2vs1$table$status <- as.factor(sign(qlf.2vs1$table$logFC))
                                            #defining names for values 1 and -1 
                                            qlf.2vs1$table$status_name <- NA
                                            #qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="1")] <- cmb_un$TARGET[[i]]
                                            #qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="-1")] <- cmb_un$SOURCE[[i]]
                                            qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="1")] <- "UP"
                                            qlf.2vs1$table$status_name[which(qlf.2vs1$table$status=="-1")] <- "DOWN"
                                            #subsetting and saving
                                            qlf.2vs1_final <- qlf.2vs1$table
                                            qlf.2vs1_final2 <- qlf.2vs1_final[which(qlf.2vs1_final$FDR < pval & abs(qlf.2vs1_final$logFC) > logFC),]
                                            
                                            message(paste0("Genes kept after filtering: ",nrow(qlf.2vs1_final2)), "printing groups...")
                                            print(tapply(qlf.2vs1_final2$status_name,qlf.2vs1_final2$status_name,length))
                                            
                                            #saving in csv tables (full and filtered)
                                            write.csv(qlf.2vs1_final, file= paste0(csv_dir,"/","glmQLFTest_",
                                                                                   cmb_un$COMP[[i]],"_","pval_",
                                                                                   as.character(pval),"_","fulltable.csv"),
                                                      row.names = T,quote = F)
                                            write.csv(qlf.2vs1_final2, file= paste0(csv_dir,"/","glmQLFTest_",
                                                                                    cmb_un$COMP[[i]],"_","pval_",
                                                                                    as.character(pval),"_","filtered.csv"),
                                                      row.names = T,quote = F) 
                                            
                                            
                                            #plot significant genes
                                            
                                            pdf(paste0(plt_dir,"/","plotMD_glmQLFTest_",cmb_un$COMP[[i]],".pdf"),height = 8, width = 8)
                                            limma::plotMD(qlf.2vs1)
                                            #abline(h=c(-1, 1), col="blue")
                                            dev.off()
                                            
                                            
                                            q1 <- data.frame(count = tapply(qlf.2vs1_final2$status_name,qlf.2vs1_final2$status_name,length))
                                            q1$groups <- row.names(q1)
                                            q1 <- q1[,c(2,1)]
                                            print(q1)
                                            #saving in csv tables (full and filtered)
                                            write.csv(q1, file= paste0(out_dir_1,"/","glmQLFTest_",
                                                                       cmb_un$COMP[[i]],"_","pval_",
                                                                       as.character(pval),"_","summary.csv"),
                                                      row.names = F)
                                            
                                            #Do not move it is need to return qlf.2vs1_final
                                            qlf.2vs1_final2
                                          })
    
    parallel::stopCluster(cl)
  } else {
    stop("no contrasts to run.") #NEW
    cont_results <- NULL # NEW
  }
  
  final_list <- list(contrasts = contrasts,
                     contrasts_results=cont_results)
  message("DONE!")
  return(final_list)
  
}

################################################################################
#parameters
#folder_name <- c(1,2,3) #folder name
pval = 0.05 # p value for filtering
logFC = 1   # logFC value for filtering 
numCores <- 4 # number of cores to use in parallel
plot_MDS <- TRUE # if plot should be appears in R session


##############################################################################
# loop over all the datasets
folder_names <- c("SRP068562", "SRP102142", "SRP132192", "SRP151473", "SRP174413", "SRP222782", "SRP062027",
                  "SRP246269", "SRP299865", "SRP338699", "SRP353076", "SRP063383", "SRP101911", "SRP125635", 
                  "SRP135093", "SRP170079", "SRP200223", "SRP226120", "SRP284189", "SRP303935", "SRP352402")
group_vect <- NULL
for(folder_name in folder_names) {
  
  ##############################################################################
  out_name <- folder_name
  ##############################################################################
  #Directories
  #where are the salmon files
  data_dir <- paste0("/scratch/grassgrns/analysis/", folder_name, "/salmon")
  #defining output folder
  out_dir <- "/ngsprojects/grassgrns/results/DEG"
  #path where the nf core metadata is available
  met_dir <-  paste0("/scratch/grassgrns/analysis/", folder_name, "/multiqc/", folder_name, "_multiqc_report_data")
  ##############################################################################
  
  x <- DEG_edgeR_func_comb(out_name,folder_name, pval,out_dir,met_dir,plot_MDS,numCores,group_vect)
  
}
