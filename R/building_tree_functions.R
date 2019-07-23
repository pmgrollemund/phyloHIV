############ ---- building_tree ---- ############
#' @title building_tree
#' @description Script for:
#' @return dqsfds
#' @param names_file_aligned_ndrm XXXX
#' @param path_output XXXX
#' @param tree_builder XXXX
#' @param bs.n XXXX
#' @param bs.id XXXX
#' @param verbose XXXX
#' @export
building_tree <- function(names_file_aligned_ndrm,path_output,
                          tree_builder="fasttree",
                          bs.n = 50,bs.id=NA,
                          verbose=TRUE){
 if(!dir.exists(path_output))
  dir.create(path_output)
 path_output <- paste(normalizePath(path_output),"/",sep="")

 if(verbose) cat("Building trees with ",tree_builder,".\n",sep="")
 for(i in 1:length(names_file_aligned_ndrm)){
  name <- names_file_aligned_ndrm[i]


  if(is.na(bs.id)){
   output_file_non_bs <- paste(path_output,tree_builder,"_",
                               gsub(x=basename(name),
                                    pattern=".fasta",
                                    replacement="/"),
                               gsub(x=basename(name),
                                    pattern=".fasta",
                                    replacement=".newick"),
                               sep="")
   if(file.exists(output_file_non_bs))
   cat('##############################################################\n',
       "The RAxML files and the output file ... \n\t",output_file_non_bs,
       "\n ... already exists, it will not be updated.\n",
       "Only the bootstrap trees will be updated.\n",
       "If you want to want to update these files, you should first remove ",
       "these files before running the current function.",
       '##############################################################\n',
       sep=""
       )
  }

  if(verbose) cat("\tFile:  ",basename(name),"\n",sep="")

  ### Make a directory for the result
  path_result <- paste(path_output,tree_builder,"_",gsub(x=basename(name),".fasta","/"),sep='')
  if(!dir.exists(path_result))
   dir.create(path_result)

  ### Name of the result file
  outfile <- paste(path_result,gsub(x=basename(name),".fasta",".newick"),sep="")

  ### Fasttree
  if(tree_builder == "fasttree"){
   if(!is.na(bs.id)){
    bs.id <- bs.id-1

    if(verbose) cat("Bootstrap id:",sprintf("%03d",bs.id),"\n")

    bs.dir <- path_result
    cmd <- cmd.fasttree.one.bootstrap(
     name,bs.id,
     outfile=file.path(bs.dir, gsub("\\.fa|\\.fasta|\\.FA|\\.FASTA",
                                    paste0("_ft.",sprintf("%03d",bs.id),".newick"),basename(name))),
     pr=big.phylo:::PR.FASTTREE, pr.args="-nt -gtr -gamma",prog.bscreate=big.phylo:::PR.EXAML.BSCREATE,
     opt.bootstrap.by="nucleotide", check.binary=TRUE)

    tmp	<- file.path(bs.dir, gsub("\\.fa|\\.fasta|\\.FA|\\.FASTA",
                                  paste0("_ft.",sprintf("%03d",0),".newick"),basename(name)))
    tmp	<- cmd.fasttree.add.bootstrap.nodelabels.to.base.tree(
     bs.dir, bs.n, tmp, outfile, bs.pattern="*_ft.[0-9][0-9][0-9].newick",
     pr=big.phylo:::PR.RAXML)

    cmd	<- paste(cmd, tmp, sep="\n")
    ### Run the commands
    system(cmd)
   }else{
    cmd <- big.phylo::cmd.fasttree.many.bootstraps(name, path_result,
                                                   bs.n, outfile,
                                                   pr.args='-nt -gtr -gamma',
                                                   opt.bootstrap.by='nucleotide')
    cmd <- strsplit(cmd,split="\n")
    n <- which(cmd[[1]] == "if [ $(find . -name '*_ft.[0-9][0-9][0-9].newick' | wc -l) -eq 2 ]; then ")
    cmd[[1]][n] <- "if [ $(find . -name '*_ft.[0-9][0-9][0-9].newick' | wc -l) -ge 2 ]; then "
    system(paste(cmd[[1]],collapse="\n"))
   }
  }

  ### ExaML
  if(tree_builder == "ExaML"){
   dname <- dirname(name)
   bname <- basename(name)

   outdir <- paste(path_tree,"ExaML_",gsub(x=bname,'.fasta',''),sep="")
   tmpdir <- paste(path_tree,"TMP_ExaML_",gsub(x=bname,'.fasta',''),sep="")
   # cmd <- paste("mkdir ",outdir," \n",sep="")
   cmd <- NULL
   cmd <- paste(cmd, unlist(big.phylo::cmd.examl.bootstrap(dname,bname,bs.to=bs.n,
                                                           outdir=outdir,
                                                           tmpdir.prefix=tmpdir)),sep="")
   cmd <- paste(cmd,"rm -R ",tmpdir,"* \n",sep="")
   cmd <- paste(cmd,"mv ",outdir,"/ExaML_result* ",outdir,"/",
                gsub(x=bname,".fasta",".newick"),"\n",sep="")
   system(cmd)
  }

  ### Maximum parsimony
  if(tree_builder == "maximum_parsimony"){
   seq <- ape::read.dna(name,format="fa")
   dist <- ape::dist.dna(seq)
   tree <- ape::nj(dist)
   ape::write.tree(tree,outfile)
  }
 }
}

############ ---- postprocess_tree ---- ############
#' @title postprocess_tree
#' @description Script for:
#' @return dqsfds
#' @param tree_names XXXX
#' @param path_outgroup XXXX
#' @param name_HXB2 XXXX
#' @param bs.n XXXX
#' @param bs.id XXXX
#' @param verbose XXXX
#' @export
postprocess_tree <- function(tree_names,
                             path_outgroup,name_HXB2,
                             bs.n = 50,bs.id=NA,
                             verbose=TRUE){
 if(verbose) cat("Postprocess the tree(s).\n ")

 for(i in 1:length(tree_names)){
  if(verbose) cat("\tFile: ",basename(tree_names[i]),"\n ",sep="")

  if(file.exists(path_outgroup)){
   load(path_outgroup)
   outgroup <- outgroup[i]
  }else{
   path_outgroup_i <-
    gsub(x=path_outgroup,
         pattern=".rda",
         replacement=paste("_",i,".rda",sep=""))
   if(file.exists(path_outgroup_i)){
    load(path_outgroup_i)
   }else{
    stop("You need to specify a path to an 'outgroup.rda' file which contains ",
         "the output of 'add_outgroup'.")
   }
  }

  if(is.na(bs.id)){
   tree_names_i <- tree_names[i]
   for(j in 1:bs.n){
    tree_names_i <-  c(
     tree_names_i,
     gsub(x= tree_names[i],
          pattern=".newick",
          replacement=paste("_ft.",sprintf("%03d",j-1),".newick",sep=""))
    )
   }
  }else{
   bs.id <- bs.id - 1
   if(file.exists(tree_names[i]) &&
      !file.exists(gsub(x=tree_names[i],".newick","_rerooted.newick"))){
    tree_names_i <- c(tree_names[i],
                      gsub(x= tree_names[i],
                           pattern=".newick",
                           replacement=paste("_ft.",sprintf("%03d",bs.id),".newick",sep=""))
    )
   }else{
    tree_names_i <- gsub(x= tree_names[i],
                           pattern=".newick",
                           replacement=paste("_ft.",sprintf("%03d",bs.id),".newick",sep="")
                         )
   }
  }

  for(name in tree_names_i){
   ### Load the tree
   tree <- ape::read.tree(name)

   ### Remove HXB2
   if(verbose) cat("\tRemoving HXB2.\n ",sep="")
   tree <- ape::drop.tip(tree,tip=which(tree$tip.label %in% name_HXB2))

   if(verbose) cat("\tRerooting.\n ",sep="")
   ### Change the tips labels
   if(sum(tree$tip.label %in% outgroup[[1]]) == 0){
    tmp2 <- strsplit(outgroup[[1]],split="[.]")[[1]][1]
    tmp <- unlist(lapply(strsplit(tree$tip.label,split="[.]"),function(v) v[1]) )

    indexes <- which(tmp == tmp2)
    tree$tip.label[indexes] <- 'OTH'
   }else{
    tree$tip.label[tree$tip.label %in% outgroup[[1]] ] <- 'OTH'
   }
   ### Cut the tree for rerooting it
   # keep only one OTH to easily cut the tree
   tmp		<- which(grepl('OTH',tree$tip.label))
   tree <- ape::drop.tip(tree,tip=tmp[-1])
   tmp  <- tmp[1]
   # and now reroot the tree
   # tree <- phytools::reroot(tree, tmp,position=median(tree$edge.length))
   tree <- ape::root(tree, tmp)
   # tree$edge.length[ which(tree$edge[,2]%in%tmp) ])
   tree <- ape::ladderize( tree )
   # and cut the last OTH
   tmp		<- which(grepl('OTH',tree$tip.label))
   # tree <- drop.tip(tree,tip=tmp)
   tree$tip.label[tmp] <- ''
   tree$edge.length[ which(tree$edge[,2] == tmp) ] <- max(tree$edge.length) / 1e10

   ape::write.tree(tree,gsub(x=name,".newick","_rerooted.newick"))
  }
 }
}
