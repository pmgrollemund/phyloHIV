############ ---- tips_relabelling ---- ############
#' @title tips_relabelling
#' @description Script for:
#' @return dqsfds
#' @param tree_names XXXX
#' @param person XXXX
#' @param sequences_meta XXXX
#' @param Country_db XXXX
#' @param vars_for_label XXXX
#' @param local_group XXX
#' @param focus_group XXXX
#' @param focus_subgroup XXX
#' @param index_subgroup XXXX
#' @param name_subgroup XXX
#' @param negative XXXX
#' @param local_indicator XXX
#' @param focus_indicator XXXX
#' @param separator XXX
#' @param blacklist_paths XXXX
#' @param blacklist_names XXX
#' @param Unknown XXXX
#' @param local_tip_split XXX
#' @param external_tip_split XXXX
#' @param local_index_location XXX
#' @param external_index_location XXXX
#' @param bs.n XXXX
#' @param bs.id XXXX
#' @param verbose XXXX
#' @export
tips_relabelling <- function(tree_names,
                             person,sequences_meta,Country_db,
                             vars_for_label,
                             local_group="WA",
                             focus_group="KC",
                             focus_subgroup=NULL,
                             index_subgroup=NULL,
                             name_subgroup=NULL,
                             negative="non",
                             local_indicator="PR/RT",
                             focus_indicator=TRUE,
                             separator="___",
                             blacklist_paths=NULL, blacklist_names=NULL,
                             Unknown = "Unknown",
                             local_tip_split="-",
                             external_tip_split="[.]",
                             local_index_location=8,
                             external_index_location=2,
                             bs.n = 50,
                             bs.id = NA,
                             verbose =TRUE){
 if(verbose) cat("Relabel the tip labels.\n")
 if(is.null(name_subgroup)) name_subgroup <- focus_subgroup

 for(i in 1:length(tree_names)){
  name_tree <- tree_names[i]

  if(is.na(bs.id)){
   names_tree_i <- name_tree
   for(j in 1:bs.n){
    names_tree_i <- c(
     names_tree_i,
     gsub(x= name_tree,
          pattern="_rerooted.newick",
          replacement=paste("_ft.",sprintf("%03d",j-1),"_rerooted.newick",sep="")
     ))
   }
  }else{
   bs.id <- bs.id - 1
   if(bs.id == 0){
    names_tree_i <- c(
     name_tree,
     gsub(x= name_tree,pattern="_rerooted.newick",
          replacement=paste("_ft.",sprintf("%03d",bs.id),"_rerooted.newick",sep="")))
   }else{
    names_tree_i <- gsub(x= name_tree,pattern="_rerooted.newick",
                         replacement=paste("_ft.",sprintf("%03d",bs.id),"_rerooted.newick",sep=""))
   }
  }

  for(name in names_tree_i){
   if(verbose) cat("\tFile:  ",basename(name),"\n",sep="")
   ### Load the tree
   tree <- ape::read.tree(name)

   ### Change some tips labels to get the same suitable format for all the tips ----
   tree <- tips_relabelling.internal(tree,
                                     person,
                                     sequences_meta,
                                     vars_for_label,
                                     local_indicator)

   ### Change the tip labels so that the group can be clearly read ----
   index <- grep(x=tree$tip.label,local_indicator)
   tmp <- tree$tip.label[index]

   # Local groups
   index2 <- which(unlist(lapply(strsplit(tmp,split=local_tip_split),
                                 function(v) v[local_index_location]) %in% focus_indicator))
   if(length(index2) > 0){
    tree$tip.label[index][index2] <-
     paste(focus_group,separator, tree$tip.label[index][index2],sep="")
   }
   if(length(index2) < length(index)){
    if(length(index2) > 0){
     tree$tip.label[index][-index2] <-
      paste(local_group,"/",focus_group,separator, tree$tip.label[index][-index2],sep="")
    }else{
     tree$tip.label[index] <-
      paste(local_group,"/",focus_group,separator, tree$tip.label[index],sep="")
    }
   }
   if(!is.null(focus_subgroup)){
    index3 <- which(unlist(lapply(strsplit(tree$tip.label[index][index2],split="-"),
                                  function(v) v[index_subgroup])) %in% focus_subgroup)
    tree$tip.label[index][index2][index3] <-
     gsub(x=tree$tip.label[index][index2][index3],
          pattern=paste(focus_group,separator,sep=""),
          replacement=paste(focus_group,name_subgroup,separator,sep=""))
    tree$tip.label[index][index2][-index3] <-
     gsub(x=tree$tip.label[index][index2][-index3],
          pattern=paste(focus_group,separator,sep=""),
          replacement=paste(focus_group,negative,name_subgroup,separator,sep=""))
   }

   # external groups
   trait <- unlist(lapply(strsplit(tree$tip.label[-index],split=external_tip_split),
                          function(v) v[external_index_location]))
   trait[trait == ""] <- Unknown
   trait[is.na(trait)] <- Unknown
   trait[trait=="-"] <- Unknown
   for(k in 1:length(Country_db)){
    trait[trait %in% Country_db[[k]]] <- names(Country_db)[k]
   }
   tree$tip.label[-index] <- paste(trait,separator, tree$tip.label[-index],sep="")

   #### The blacklisted: node index for those the state is unknown ----
   blacklist <- tree$tip.label[grep(x=tree$tip.label,paste(Unknown,separator,"[.]",sep=""))]

   #### Save the results ----
   # The tree with the new labels
   new_name <- gsub(x=name,pattern=".newick",replacement="_renamed.newick")
   ape::write.tree(tree,new_name)

   # The Blacklist file
   if(is.null(blacklist_paths)){
    blacklist_path <- dirname(name)
   }else{
    blacklist_path <- blacklist_paths[i]
   }
   if(is.null(blacklist_names)){
    blacklist_name <- paste("Blacklist_",
                            gsub(x=basename(name),
                                 pattern=".newick",
                                 replacement="_renamed.txt"
                            ),
                            sep="")
   }else{
    blacklist_name <- blacklist_names[i]
   }
   write(blacklist,paste(blacklist_path,"/",blacklist_name,sep=""))
  }
 }
}

############ ---- tips_relabelling.internal ---- ############
#' @title tips_relabelling.internal
#' @description Script for:
#' @return dqsfds
#' @param tree XXXX
#' @param person XXXX
#' @param sequences_meta XXX
#' @param vars_for_label XXXX
#' @param local_indicator XXXX
#' @export
tips_relabelling.internal <- function(tree,person,sequences_meta,vars_for_label,
                                      local_indicator="PR/RT"){
 indexes <- grep(x=tree$tip.label,local_indicator)

 meta_tmp <- sequences_meta[sequences_meta$seqID %in% tree$tip.label[indexes],
                            c("newnum","seqID")]
 IDS      <- as.numeric(meta_tmp$newnum)

 index_newnum <- which(names(person) == "newnum")
 index_vars <- which(names(person) %in% vars_for_label)

 labels <- data.table::data.table("newnum" = person[person$newnum %in% IDS,index_newnum],
                      "label"  = apply(person[person$newnum %in% IDS,index_vars],1,
                                     function(v) paste(v,collapse="-")))
 names(labels) <- c("newnum","label")

 # data.table::setkey(meta_tmp,newnum)
 # data.table::setkey(labels,newnum)
 # labels <- labels[meta_tmp]
 # The previous commands (which depend on the 'data.table' package) do not
 # work with the package.
 # The following commands replace it.
 # Begin:
 order_newnum <- sapply(labels$newnum, function(n) which(meta_tmp$newnum == n))
 meta_tmp <- meta_tmp[order_newnum,]
 labels$seqID <- meta_tmp$seqID
 # End.

 labels$label <- gsub(x=labels$label," ","")
 labels$label <- paste(labels$seqID,labels$label,sep="-")
 labels$label <- gsub(x=labels$label,"--","-_-")
 data.table::setkey(labels,seqID)

 order_tmp <- order(tree$tip.label[indexes])

 labels$label <- gsub(x=labels$label,"'","")
 tree$tip.label[indexes][order_tmp] <- labels$label

 return(tree)
}

############ ---- wrap_phyloscanner_analyse_tree ---- ############
#' @title wrap_phyloscanner_analyse_tree
#' @description Script for:
#' @return dqsfds
#' @param tree_names XXXX
#' @param path_output XXXX
#' @param opt_regex XXXX
#' @param bs.n XXXX
#' @param bs.id XXXX
#' @param verbose XXX
#' @export
wrap_phyloscanner_analyse_tree <- function(tree_names,path_output,
                                           opt_regex="([[:alpha:]]+)___.*",
                                           bs.n = 50,bs.id = NA,
                                           verbose =TRUE){
 path_phyloscanner <- system.file(package="phyloHIV","ext","phyloscanner_analyse_trees.R")
 if(verbose) cat("Estimate the ancestral states.\n")
 if(!dir.exists(path_output))
  dir.create(path_output)

 for(i in 1:length(tree_names)){
  name_tree <- tree_names[i]
  if(verbose) cat("\tFile:  ",basename(name_tree),"\n",sep="")

  path_res <- paste(path_output,gsub(x=basename(name_tree),".newick",""),sep="")
  if(!dir.exists(path_res))
   dir.create(path_res)

  if(is.na(bs.id)){
   names_tree_i <- name_tree
   for(j in 1:bs.n){
    names_tree_i <- c(
     names_tree_i,
     gsub(x= name_tree,
          pattern="_rerooted_renamed.newick",
          replacement=paste("_ft.",sprintf("%03d",j-1),"_rerooted_renamed.newick",sep="")
     ))
   }
  }else{
   bs.id <- bs.id - 1
   if(bs.id == 0){
    names_tree_i <- c(
     name_tree,
     gsub(x= name_tree,pattern="_rerooted_renamed.newick",
          replacement=paste("_ft.",sprintf("%03d",bs.id),"_rerooted_renamed.newick",sep="")))
   }else{
    names_tree_i <- gsub(x= name_tree,pattern="_rerooted_renamed.newick",
                         replacement=paste("_ft.",sprintf("%03d",bs.id),"_rerooted_renamed.newick",sep=""))
   }
  }

  for(name in names_tree_i){
   blacklist_path <- paste(dirname(name),"/Blacklist_",
                           gsub(x=basename(name),
                                pattern=".newick",
                                replacement=".txt"),sep="")

   options <- NULL
   options <- paste(options,name,sep="") # input path
   options <- paste(options," ",basename(name),sep="") # job title
   options <- paste(options,' s,20 -m 1e-5 -x "',opt_regex,'" -v 2 -ow -rda ',sep="")
   options <- paste(options,' -od ',path_res,sep="") # output plot name
   options <- paste(options," -b ",blacklist_path,sep="") # path to blacklisted nodes index
   cmd <- paste("Rscript",path_phyloscanner,options,sep=" ")

   # Run the command
   system(cmd)
  }
 }
}

############ ---- phyloscanner_postprocess ---- ############
#' @title phyloscanner_postprocess
#' @description Script for:
#' @return dqsfds
#' @param tree_names XXXX
#' @param path_output XXXX
#' @param focus_group XXXX
#' @param bs.n XXXX
#' @param bs.id XXXX
#' @param verbose XXX
#' @export
phyloscanner_postprocess <- function(tree_names,path_output,
                                     focus_group,
                                     bs.n=50,bs.id=NA,
                                     verbose=TRUE){
 if(verbose) cat("Postprocess the colored trees.\n")
 for(i in 1:length(tree_names)){
  name_tree <- tree_names[i]
  if(verbose) cat("\tFile:  ",basename(name_tree),"\n",sep="")

  path_res <- paste(path_output,gsub(x=basename(name_tree),".newick",""),sep="")

  if(is.na(bs.id)){
   names_tree_i <- name_tree
   for(j in 1:bs.n){
    names_tree_i <- c(
     names_tree_i,
     gsub(x= name_tree,
          pattern="_rerooted_renamed.newick",
          replacement=paste("_ft.",sprintf("%03d",j-1),"_rerooted_renamed.newick",sep="")
     ))
   }
  }else{
   bs.id <- bs.id - 1
   if(bs.id == 0){
    names_tree_i <- c(
     name_tree,
     gsub(x= name_tree,pattern="_rerooted_renamed.newick",
          replacement=paste("_ft.",sprintf("%03d",bs.id),"_rerooted_renamed.newick",sep="")))
   }else{
    names_tree_i <- gsub(x= name_tree,pattern="_rerooted_renamed.newick",
                         replacement=paste("_ft.",sprintf("%03d",bs.id),"_rerooted_renamed.newick",sep=""))
   }
  }

  for(name in names_tree_i){
   res_name <- paste(basename(name),"_workspace.rda",sep="")
   load(paste(path_res,"/",res_name,sep=""))

   tree <- phyloscanner.trees[[1]]$tree
   tmp <- phyloscanner.trees[[1]]$clade.mrcas.by.host
   index <- which(names(tmp) == focus_group)
   if(length(index)==1){
    MRCA <- tmp[[index]]
   }else{
    MRCA <- NULL
   }

   # Some fixes
   tmp <- as.character(attributes(tree)$BRANCH_COLOURS)
   tmp[is.na(tmp)] <- 'Unknown'
   index_focus_group <- grep(x=tree$tip.label,paste("^",focus_group,sep=""))
   tmp[index_focus_group] <- focus_group
   trait <- data.frame(BRANCH_COLOURS = tmp,
                       INDIVIDUAL     = as.character(attributes(tree)$INDIVIDUAL),
                       SUBGRAPH_MRCA  = as.character(attributes(tree)$SUBGRAPH_MRCA))

   # Save the files
   # write.tree(tree,paste(path_tree,
   #                       "fasttree_newseq_LANL_Subtype_",subtype,"_mafft_aligned_ndrm/",
   #                       "newseq_LANL_Subtype_",subtype,"_mafft_aligned_ndrm_rerooted_renamed_mp.newick",sep=""))
   utils::write.csv(trait,paste(path_res,"/",
                         gsub(x=basename(name),
                              pattern=".newick",
                              replacement="_trait.csv"),
                         sep=""),
             row.names=F)
   utils::write.csv(as.character(MRCA),paste(path_res,"/",
                                      gsub(x=basename(name),
                                           pattern=".newick",
                                           replacement="_MRCA.txt"),
                                      sep=""),
             row.names=F)
  }
 }
}

############ ---- clades_extracting ---- ############
#' @title clades_extracting
#' @description Script for:
#' @return dqsfds
#' @param tree_names XXXX
#' @param path_colored_tree XXXX
#' @param path_output XXXX
#' @param focus_group XXXX
#' @param bs.n XXXX
#' @param bs.id XXXX
#' @param verbose XXX
#' @export
clades_extracting <- function(tree_names,
                              path_colored_tree,path_output,
                              focus_group,
                              bs.n=50,bs.id=NA,
                              verbose=T){
 if(verbose) cat("Extract the clades.\n")

 if(!dir.exists(path_output))
  dir.create(path_output)

 if(is.na(bs.id)){
  bs.id_index <- 0:bs.n
 }else{
  bs.id_index <- bs.id
  if(bs.id_index==1)
   bs.id_index <- c(0,bs.id_index)
 }

 for(bs.id in bs.id_index){
  clades <- list()

  for(i in 1:length(tree_names)){
   name_tree <- tree_names[i]

   if(bs.id!=0){
    name_tree <- gsub(x= name_tree,
                      pattern="_rerooted_renamed.newick",
                      replacement=paste("_ft.",sprintf("%03d",bs.id-1),"_rerooted_renamed.newick",sep=""))
   }

   name_MRCA_file <- paste(path_colored_tree,
                           gsub(x=basename(tree_names[i]),
                                pattern=".newick",
                                replacement=""),
                           "/",
                           gsub(x=basename(name_tree),
                                pattern=".newick",
                                replacement="_MRCA.txt"),
                           sep="")

   if(!file.exists(name_MRCA_file)){
    cat("All the bootstrap are not computed. Wait for it... \n")
    return()
   }


   if(verbose) cat("\tFile:  ",basename(name_tree),"\n",sep="")


   #### Load the tree, trait and MRCA
   tree <- ape::read.tree(name_tree)
   trait <- utils::read.csv(paste(path_colored_tree,
                                  gsub(x=basename(tree_names[i]),
                                       pattern=".newick",
                                       replacement=""),
                                  "/",
                                  gsub(x=basename(name_tree),
                                       pattern=".newick",
                                       replacement="_trait.csv"),
                                  sep=""),header=T)

   if(file.size(name_MRCA_file) > 1){
    MRCA <- utils::read.csv(name_MRCA_file,header=T)[,1]
    MRCA <- MRCA[!is.na(MRCA)]

    if(length(MRCA)>0){
     ### MRCA correction
     MRCA <- MRCA_correction(tree,trait,MRCA,focus_group)

     #### Clades extracting
     tmp_clades <- wrap_extract.clades(tree,trait,MRCA) # see the function in XXXX
     names(tmp_clades) <- rep(basename(name_tree),length(tmp_clades)) # The names of the list items are the subtypes
     clades <- c(clades,tmp_clades)
    }
   }
  }

  if(bs.id==0){
   save(clades,file=paste(path_output,"Clades_save.RData",sep=""))
  }else{
   save(clades,file=paste(path_output,"Clades_save_",sprintf("%03d",bs.id-1),".RData",sep=""))
  }

 }
}
############ ---- MRCA_correction ---- ############
#' @title clades_extracting
#' @description Script for:
#' @return dqsfds
#' @param tree XXXX
#' @param traits XXXX
#' @param MRCAs XXXX
#' @param state XXXX
#' @export
MRCA_correction <- function(tree,traits,MRCAs,state){
 tree$node.label <- traits$BRANCH_COLOURS[1:length(tree$node.label) + length(tree$tip.label)]

 anc_MRCA <- phangorn::Ancestors(tree,MRCAs)
 sub_clades <- which(traits[unlist(lapply(anc_MRCA,function(v) v[1])),1] == state )

 if(length(sub_clades) > 0){
  new_MRCA <- NULL
  for(i in 1:length(sub_clades)){
   anc_traits <- traits[anc_MRCA[[sub_clades[i]]],1]
   anc_index  <- rev(which(anc_traits == state))[1]

   new_MRCA <-  c(new_MRCA,anc_MRCA[[sub_clades[i]]][anc_index])
  }
  MRCAs <- c(MRCAs[-sub_clades],unique(new_MRCA))
 }

 return(MRCAs)
}

############ ---- extract.clades_wrap ---- ############
#' @title extract.clades_wrap
#' @description Script for:
#' @return dqsfds
#' @param tree XXXX
#' @param traits XXXX
#' @param MRCAs XXXX
#' @export
wrap_extract.clades <- function(tree,traits,MRCAs)
{
 tree$node.label <- traits$BRANCH_COLOURS[1:length(tree$node.label) + length(tree$tip.label)]

 clades <- list()
 Unselected_tips <- NULL
 for(i in 1:length(MRCAs)){

  if(MRCAs[i] < (length(tree$tip.label) + tree$Nnode)){
   if(MRCAs[i] > length(tree$tip.label)){
    clade_tmp <- ape::extract.clade(tree,MRCAs[i])

    clade_anc <- phangorn::Ancestors(tree,MRCAs[i])[1]
    clade_tmp$source_location <- tree$node.label[clade_anc-length(tree$tip.label)]
    clade_tmp$source_location[is.na(clade_tmp$source_location)] <- "Unknown"
    if(length(clade_tmp$source_location) > 0){
     if(clade_tmp$source_location == "Unknown"){
      trait_ancs <- tree$node.label[phangorn::Ancestors(tree,MRCAs[i])-length(tree$tip.label)]
      first_anc_trait <- which(trait_ancs != 'Unknown')[1]
      if(!is.na(first_anc_trait))
       clade_tmp$source_location <- trait_ancs[first_anc_trait]
     }
    }
    clade_tmp$length_root <- tree$edge.length[clade_anc]

    clades <- c(clades,list(clade_tmp))
   }else{
    Unselected_tips <- c(Unselected_tips,MRCAs[i])
   }

  }
 }

 for(i in Unselected_tips){
  other_tips <- (1:length(tree$tip.label))[-i]
  clade_tmp <- ape::drop.tip(tree,tip=other_tips)

  clade_anc <- phangorn::Ancestors(tree,i)[1]
  clade_tmp$source_location <- tree$node.label[clade_anc-length(tree$tip.label)]
  clade_tmp$source_location[is.na(clade_tmp$source_location)] <- "Unknown"

  if(length(clade_tmp$source_location) > 0){
   if(clade_tmp$source_location == "Unknown"){
    trait_ancs <- tree$node.label[phangorn::Ancestors(tree,i)-length(tree$tip.label)]
    first_anc_trait <- which(trait_ancs != 'Unknown')[1]
    if(!is.na(first_anc_trait))
     clade_tmp$source_location <- trait_ancs[first_anc_trait]
   }
  }
  clade_tmp$length_root <- tree$edge.length[clade_anc]

  clades <- c(clades,list(clade_tmp))
 }

 return(clades)
}
