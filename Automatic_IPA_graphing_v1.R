require(ggplot2)
require(reshape2)
require(stringr)

## Description:
## An R script for generating summary plots of parsed IPA files generated from
## "Parse_IPA.py". At the moment, this R script only plots the results from
## three of the IPA result types. Those three are "Upstream regulators", 
## "Causal Networks", and "Canonical Pathways".
## Should a user wish to adapt the script to also generate graphs for 
## "Diseases and Bio Functions" and "Tox Functions", it should be simple enough to
## adapt the current script for the additional graphs. 

## Helper functions
ReadFileForIPAplot <- function(filename_prefix){
  data_upreg <- read.delim(paste0(filename_prefix,"_upstreamregulators.txt"), sep="\t", stringsAsFactors = F, header = T)
  data_causal <- read.delim(paste0(filename_prefix,"_causalnetworks.txt"), sep="\t", stringsAsFactors = F, header = T)
  data_canon <- read.delim(paste0(filename_prefix,"_canonicalpathways.txt"), sep="\t", stringsAsFactors = F, header = T)
  # try(data_dis <- read.delim(paste0(filename_prefix,"_diseases.txt"), sep="\t", stringsAsFactors = F, header = T),silent=TRUE)
  # try(data_tox <- read.delim(paste0(filename_prefix,"_toxfunctions.txt"), sep="\t", stringsAsFactors = F, header = T),silent=TRUE)
  return(list(data_upreg,data_causal,data_canon))
}

fix_slots <- function(x,curmat,
                      select_by="Upstream.Regulator",
                      extract_by="Bias.corrected.z.score"){
  intermat <- curmat[curmat[,select_by]==x,]
  if(all(is.na(intermat[,extract_by]))){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

## Supplement missing entries
clean_df <- function(curdf,feat_vec,select_by="Upstream.Regulator",
                     cat_1="Analysis",fix_by="Bias.corrected.z.score",
                     replace_vec=c(0),current_set=c_set){
  out_df <- curdf
  # loop through features
  for(curvar in feat_vec){
    d5 <- curdf[curdf[,select_by]==curvar,]
    if(all(current_set%in%d5[,cat_1])){
      out_df <- rbind(d5,out_df)
    }else{
      ## which features are missing?
      missing_feat <- current_set[!(current_set%in%d5[,cat_1])]
      for(cur_miss in missing_feat){
        d5[nrow(d5)+1,c(cat_1,select_by,fix_by)] <- c(cur_miss,curvar,replace_vec)
      }
      out_df <- rbind(d5,out_df)
    }
  }
  ## remove duplicate entries
  out_df <- out_df[!(duplicated(out_df[,c(cat_1,select_by)])),]
  return(out_df)
}
## End helper functions


## Primary function for generating plots
gen_plot <- function(c_set,k=3,canon_select=NULL,
                     causal_select=NULL,upreg_select=NULL,
                     pdf_prefix="IPA_plots_top3"){
  ## designate the input comparisons
  c_list=lapply(c_set,function(x) b[x])
  ## specify the output figure file
  out_file=paste0(pdf_prefix,".pdf")
  print(paste0("Generating figure: ",out_file))
  pdf(out_file)
  ## compile the upstream regulators
  d=lapply(c_list,function(x) x[[1]][[1]])
  for(i in 1:length(c_set)){
    d[[i]]$Analysis=c_set[i]
    ## Remove NA entries
    d[[i]]=d[[i]][!(is.na(d[[i]]$p.value.of.overlap)),]
  }
  d2=do.call(rbind,d)
  ## Identify the top k enriched features
  for(i in 1:length(c_set)){
    d3=d[[i]]
    d3=d3[order(d3$p.value.of.overlap,decreasing=FALSE),] ## sort the dataframe
    if(i==1){k_names=d3[1:k,"Upstream.Regulator"]}else{k_names=c(k_names,d3[1:k,"Upstream.Regulator"])}
    if(i==1){sig_tot=nrow(d3[d3$p.value.of.overlap<0.01,])}else{
      sig_tot=c(sig_tot,nrow(d3[d3$p.value.of.overlap<0.01,]))}
  }
  k_names <- sort(unique(k_names))
  names(sig_tot) <- c_set
  ##Filter out features with Z-score==NA in all entries
  if(is.null(causal_select)){
    d3 <- d2[d2$Upstream.Regulator%in%k_names,]
    z_valid <- sapply(k_names,fix_slots,d3,"Upstream.Regulator","Bias.corrected.z.score")
  }else{
    d3 <- d2[d2$Upstream.Regulator%in%causal_select,]
    z_valid <- sapply(causal_select,fix_slots,d3,select_by="Upstream.Regulator",extract_by="Bias.corrected.z.score")
  }
  d4 <- d3[d3$Upstream.Regulator%in%names(z_valid[z_valid]),]
  ## Supplement missing entries
  d3 <- clean_df(d3,feat_vec=unique(d3$Upstream.Regulator),
                 select_by="Upstream.Regulator",
                 cat_1="Analysis",
                 fix_by=c("Bias.corrected.z.score","p.value.of.overlap"),
                 replace_vec=c(0,1),current_set=c_set)
  d3$Bias.corrected.z.score <- as.numeric(d3$Bias.corrected.z.score)
  d3$p.value.of.overlap <- as.numeric(d3$p.value.of.overlap)
  
  d4 <- clean_df(d4,feat_vec=unique(d3$Upstream.Regulator),
                 select_by="Upstream.Regulator",
                 cat_1="Analysis",
                 fix_by=c("Bias.corrected.z.score","p.value.of.overlap"),
                 replace_vec=c(0,1),current_set=c_set)
  d4$Bias.corrected.z.score <- as.numeric(d4$Bias.corrected.z.score)
  d4$p.value.of.overlap <- as.numeric(d4$p.value.of.overlap)
  d3$Analysis <- factor(d3$Analysis,levels=c_set)
  d4$Analysis <- factor(d4$Analysis,levels=c_set)
  
  p <- ggplot(d4,aes(x=Bias.corrected.z.score,y=Upstream.Regulator,fill=Analysis))+
    geom_bar(stat='identity',position='dodge',width=0.7)+theme_bw()+
    theme(axis.text=element_text(color='black'))
  print(p)
  
  p <- ggplot(d3,aes(x=-log10(p.value.of.overlap),y=Upstream.Regulator,fill=Analysis))+
    geom_bar(stat='identity',position='dodge',width=0.7)+theme_bw()+
    theme(axis.text=element_text(color='black'))
  print(p)
  ## Plot the total # of significant entries
  sig_tot_mt <- melt(sig_tot)
  sig_tot_mt$id <- rownames(sig_tot_mt)
  p <- ggplot(sig_tot_mt,aes(x=id,y=value))+geom_bar(stat='identity')+
    theme_bw()+theme(axis.text=element_text(color='black'),
                     axis.text.x=element_text(angle=45,hjust=1))+
    xlab("Datasets")+ylab("Significantly enriched Upstream Regulators")
  print(p)
  
  
  
  ## compile the causal networks  (Master.Regulator) Activation.z.score
  d=lapply(c_list,function(x) x[[1]][[2]])
  for(i in 1:length(c_set)){
    d[[i]]$Analysis=c_set[i]
    ## Remove NA entries
    d[[i]]=d[[i]][!(is.na(d[[i]]$p.value.of.overlap)),]
  }
  d2=do.call(rbind,d)
  ## Identify the top k enriched features
  # k_names=c()
  for(i in 1:length(c_set)){
    d3=d[[i]]
    d3=d3[order(d3$p.value.of.overlap,decreasing=FALSE),] ## sort the dataframe
    if(i==1){k_names=d3[1:k,"Master.Regulator"]}else{k_names=c(k_names,d3[1:k,"Master.Regulator"])}
    if(i==1){sig_tot=nrow(d3[d3$p.value.of.overlap<0.01,])}else{
      sig_tot=c(sig_tot,nrow(d3[d3$p.value.of.overlap<0.01,]))}
  }
  k_names <- sort(unique(k_names))
  names(sig_tot) <- c_set
  if(is.null(causal_select)){
    d3 <- d2[d2$Master.Regulator%in%k_names,]
    z_valid <- sapply(k_names,fix_slots,d3,"Master.Regulator","Activation.z.score")
  }else{
    d3 <- d2[d2$Master.Regulator%in%causal_select,]
    z_valid <- sapply(causal_select,fix_slots,d3,"Master.Regulator","Activation.z.score")
  }
  ##Filter out features with Z-score==NA in all entries
  d4 <- d3[d3$Master.Regulator%in%names(z_valid[z_valid]),]
  
  ## Clean features
  d3 <- clean_df(d3,feat_vec=unique(d3$Master.Regulator),
                       select_by="Master.Regulator",
                       cat_1="Analysis",
                       fix_by=c("Activation.z.score","p.value.of.overlap"),
                       replace_vec=c(0,1),current_set=c_set)
  d3$Activation.z.score <- as.numeric(d3$Activation.z.score)
  d3$p.value.of.overlap <- as.numeric(d3$p.value.of.overlap)
  
  d4 <- clean_df(d4,feat_vec=unique(d3$Master.Regulator),
                 select_by="Master.Regulator",
                 cat_1="Analysis",
                 fix_by=c("Activation.z.score","p.value.of.overlap"),
                 replace_vec=c(0,1),current_set=c_set)
  d4$Activation.z.score <- as.numeric(d4$Activation.z.score)
  d4$p.value.of.overlap <- as.numeric(d4$p.value.of.overlap)
  
  d3$Analysis <- factor(d3$Analysis,levels=c_set)
  d4$Analysis <- factor(d4$Analysis,levels=c_set)
  
  ## Plot the Z-score for causal networks
  p <- ggplot(d4,aes(x=Activation.z.score,y=Master.Regulator,fill=Analysis))+
    geom_bar(stat='identity',position='dodge',width=0.7)+theme_bw()+
    theme(axis.text=element_text(color='black'))
  print(p)
  ## Plot the P-value for causal networks
  p <- ggplot(d3,aes(x=-log10(p.value.of.overlap),y=Master.Regulator,fill=Analysis))+
    geom_bar(stat='identity',position='dodge',width=0.7)+theme_bw()+
    theme(axis.text=element_text(color='black'))
  print(p)
  ## Plot the total # of significant entries for causal networks
  sig_tot_mt <- melt(sig_tot)
  sig_tot_mt$id <- rownames(sig_tot_mt)
  p <- ggplot(sig_tot_mt,aes(x=id,y=value))+geom_bar(stat='identity')+
    theme_bw()+theme(axis.text=element_text(color='black'),
                     axis.text.x=element_text(angle=45,hjust=1))+
    xlab("Datasets")+ylab("Significantly enriched Master Regulators")
  print(p)
  
  ## compile the canonical pathways
  d=lapply(c_list,function(x) x[[1]][[3]])
  for(i in 1:length(c_set)){
    d[[i]]$Analysis=c_set[i]
    ## Remove NA entries
    d[[i]]=d[[i]][!(is.na(d[[i]]$X.log.p.value.)),]
  }
  d2=do.call(rbind,d)
  ## Identify the top k enriched features
  for(i in 1:length(c_set)){
    d3=d[[i]]
    d3=d3[order(d3$X.log.p.value.,decreasing=TRUE),] ## sort the dataframe
    if(i==1){k_names=d3[1:k,"Ingenuity.Canonical.Pathways"]}else{k_names=c(k_names,d3[1:k,"Ingenuity.Canonical.Pathways"])}
    if(i==1){sig_tot=nrow(d3[d3$X.log.p.value.<0.01,])}else{
      sig_tot=c(sig_tot,nrow(d3[d3$X.log.p.value.<0.01,]))}
  }
  k_names <- sort(unique(k_names))
  names(sig_tot) <- c_set
  ##Filter out features with Z-score==NA in all entries
  if(is.null(canon_select)){
    d3 <- d2[d2$Ingenuity.Canonical.Pathways%in%k_names,]
    z_valid <- sapply(k_names,fix_slots,d3,"Ingenuity.Canonical.Pathways","zScore")
  }else{
    d3 <- d2[d2$Ingenuity.Canonical.Pathways%in%canon_select,]
    z_valid <- sapply(canon_select,fix_slots,d3,"Ingenuity.Canonical.Pathways","zScore")
  }
  d4 <- d3[d3$Ingenuity.Canonical.Pathways%in%names(z_valid[z_valid]),]
  ## Clean features
  d3 <- clean_df(d3,feat_vec=unique(d3$Ingenuity.Canonical.Pathways),
                 select_by="Ingenuity.Canonical.Pathways",
                 cat_1="Analysis",
                 fix_by=c("zScore","X.log.p.value."),
                 replace_vec=c(0,1),current_set=c_set)
  d3$zScore <- as.numeric(d3$zScore)
  d3$X.log.p.value. <- as.numeric(d3$X.log.p.value.)
  
  d4 <- clean_df(d4,feat_vec=unique(d3$Ingenuity.Canonical.Pathways),
                 select_by="Ingenuity.Canonical.Pathways",
                 cat_1="Analysis",
                 fix_by=c("zScore","X.log.p.value."),
                 replace_vec=c(0,1),current_set=c_set)
  d4$zScore <- as.numeric(d4$zScore)
  d4$X.log.p.value. <- as.numeric(d4$X.log.p.value.)  

  p <- ggplot(d4,aes(x=zScore,y=Ingenuity.Canonical.Pathways,fill=Analysis))+
    geom_bar(stat='identity',position='dodge',width=0.7)+theme_bw()+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40))+
    theme(axis.text=element_text(color='black'))
  print(p)
  
  p <- ggplot(d3,aes(x=X.log.p.value.,y=Ingenuity.Canonical.Pathways,fill=Analysis))+
    geom_bar(stat='identity',position='dodge',width=0.7)+theme_bw()+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40))+
    theme(axis.text=element_text(color='black'))
  print(p)
  ## Plot the total # of significant entries
  sig_tot_mt <- melt(sig_tot)
  sig_tot_mt$id <- rownames(sig_tot_mt)
  p <- ggplot(sig_tot_mt,aes(x=id,y=value))+geom_bar(stat='identity')+
    theme_bw()+theme(axis.text=element_text(color='black'),
                     axis.text.x=element_text(angle=45,hjust=1))+
    xlab("Datasets")+ylab("Significantly enriched Canonical Pathways")
  print(p)
  
  dev.off()
}
## End main function


##Designate output directory
out_dir="C:\\Users\\Molly\\Documents\\ipa"
setwd(out_dir)

## Specify your input files
a=list("a.vs.ctl_ko"="C:\\Users\\Molly\\Documents\\ipa\\ipa_treatA.vs.control_knockout.txt",
       "a.vs.ctl_wt"="C:\\Users\\Molly\\Documents\\ipa\\ipa_treatA.vs.control_wildtype.txt",
       "b.vs.ctl_ko"="C:\\Users\\Molly\\Documents\\ipa\\ipa_treatB.vs.control_knockout.txt",
       "b.vs.ctl_wt"="C:\\Users\\Molly\\Documents\\ipa\\ipa_treatB.vs.control_wildtype.txt")
b=lapply(a,function(x) ReadFileForIPAplot(x[1]))

## Optional - Example IPA features for demonstration purposes
# upreg_example=c("l-asparaginase","CEBPB","TP53","E2F4","Eldr","CSF2",
#                 "CDKN1A","torkinib","CKAP2L","lipopolysaccharide",
#                 "dextran sulfate","MYC")
# causal_example=c("GS-444217","CEBPB","l-asparaginase","ETS1","LMO2","Eldr",
#                  "DACH1","SMARCB1","CTDSP1","1-octanol","RPL11","CYP2B6",
#                  "FBXL14","MLXIPL","SAFB","SAFB2" )
# canon_example=c("Acetone Degradation I (to Methylglyoxal)",                    
#                 "Activation of IRF by Cytosolic Pattern Recognition Receptors",
#                 "Agrin Interactions at Neuromuscular Junction",                
#                 "Altered T Cell and B Cell Signaling in Rheumatoid Arthritis", 
#                 "B Cell Development","B Cell Receptor Signaling",                                   
#                 "BAG2 Signaling Pathway","CDX Gastrointestinal Cancer Signaling Pathway",               
#                 "CNTF Signaling","CREB Signaling in Neurons")       

## Specify which IPA results you want to compare
ipa_comparisons_to_make <- c("a.vs.ctl_ko","b.vs.ctl_ko")
gen_plot(ipa_comparisons_to_make,pdf_prefix="Plot_1")

ipa_comparisons_to_make <- c("a.vs.ctl_ko","a.vs.ctl_wt",
                             "b.vs.ctl_ko","b.vs.ctl_wt")
gen_plot(ipa_comparisons_to_make,k=5,pdf_prefix="Plot_2")

