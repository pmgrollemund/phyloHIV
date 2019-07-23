############ ---- makeblastdb ---- ############
#' @title makeblastdb
#' @description Script for:
#' @return dqsfds
#' @param LANL_seq XXXX
#' @param output_file XXXX
#' @param index XXX
#' @param verbose XXXX
#' @export
makeblastdb <- function(LANL_seq,output_file,
                        index=NA,
                        verbose=T
){
 if(index == "NA") index <- NA

 ### Stop if wrong args
 stopifnot( length(LANL_seq) == length(output_file) )

 ### Which files ?
 if(is.na(index)){
  index <- 1:length(LANL_seq)
 }

 ### For each item
 if(verbose) cat("Make a Blast database.\n")
 for(i in index){
  if(i > length(LANL_seq)){
   warning("The index job is too large.\n")
  }else{
   if(verbose) cat(paste("\t",LANL_seq[i],".\n",sep=""))
   cmd <- paste("makeblastdb -in ",LANL_seq[i]," ",
                "-dbtype nucl -title ",gsub(".fasta","",basename(LANL_seq[i]))," ",
                "-out ",output_file[i],sep="")

   system(command=cmd,ignore.stdout=T)
  }
 }
}
############ ---- blastn ---- ############
#' @title blastn
#' @description Script for:
#' @return dqsfds
#' @param names_seq XXXX
#' @param names_db XXXX
#' @param nbre_close XXXX
#' @param verbose XXXX
#' @param index XXX
#' @export
blastn <- function(names_seq,names_db,
                   nbre_close=20,
                   index=NA,
                   verbose=TRUE
){
 if(index == "NA") index <- NA

 ### Stop if wrong args
 stopifnot( length(names_seq) == length(names_db) )

 ### Which files ?
 if(is.na(index)){
  index <- 1:length(names_seq)
 }

 ### For each item
 if(verbose) cat("Run blastn.\n")
 for(i in index){
  if(i > length(names_seq)){
   warning("The index job is too large.\n")
  }else{
   if(verbose) cat(paste("\t",names_seq[i],".\n",sep=""))
   cmd <- paste("blastn -query ",names_seq[i]," -db ",names_db[i], " -out ",
                gsub(".db","_closest.txt",names_db[i]),
                " -max_target_seqs ",nbre_close," -outfmt 6",sep="")
   system(command=cmd,ignore.stdout=F,ignore.stderr=F )
  }
 }
}
############ ---- get_closest_seq ---- ############
#' @title get_closest_seq
#' @description Script for:
#' @return dqsfds
#' @param LANL_seq XXXX
#' @param Local_seq XXXX
#' @param path_output XXXX
#' @param output_names XXXX
#' @param max_closest XXXX
#' @param min_percent XXXX
#' @param nbre_close XXXX
#' @param name_HXB2_file XXXX
#' @param index XXXX
#' @param verbose XXXX
#' @export
get_closest_seq <- function(LANL_seq,Local_seq,
                            path_output,output_names,
                            max_closest,
                            nbre_close=NULL,
                            min_percent,
                            name_HXB2_file=NULL,
                            index=NA,
                            verbose=TRUE){
 if(index == "NA") index <- NA

 # LANL_seq <- gsub(x=LANL_seq,
 #                  pattern="_stripped.fasta",
 #                  replacement=".fasta")

 if(verbose) cat("Get the closest LANL sequences and merge them with the local sequences.\n")
 prop <- NULL
 nbre_close_res <- NULL

 if(is.null(path_output)){
  path_output <- c(dirname(LANL_seq),"../Closest/",sep="")
 }
 if(!dir.exists(path_output)){
  dir.create(path_output)
 }

 ### Which files ?
 if(is.na(index)){
  index <- 1:length(LANL_seq)
 }

 for(i in index){
  if(i > length(LANL_seq))
   stop("The index job is too large.\n")

  name <- LANL_seq[i]
  if(verbose) cat("\tFile: ",basename(name),"\n",sep="")
  ### Initialize nbre_close_tmp
  if(is.null(nbre_close)) nbre_close_tmp <- max_closest else nbre_close_tmp <- nbre_close

  ### Load LANL data and closest file
  LANL_seq_dna <- ape::read.dna(name,format="fa")
  closest  <- data.table::as.data.table(
   utils::read.csv(gsub(x=name,pattern=".fasta","_closest.txt"),
                   header=F,sep="\t",stringsAsFactors=FALSE))
  data.table::setnames(closest, paste('V',1:12,sep=''),
                       c('new_seq_ID','LANL_ID','pident','length','mismatch','gapopen',
                         'qstart','qend','sstart','send','evalue','bitscore'))

  ### Choose the 'nbre_close' closest sequences
  if(!is.null(nbre_close)){
   if(verbose) cat("\tKeep only the ",nbre_close," closest sequences.\n",sep="")
   index <- rep(1:nbre_close,length(unique(closest$new_seq_ID))) +
    rep(1:length(unique(closest$new_seq_ID))-1,each=nbre_close) * max_closest
   # cat("bulu 1 \n")
   closest<- closest[index, ]
  }

  ### Determine the closest
  match        <- base::subset(closest,select=c(new_seq_ID,LANL_ID))
  closest_LANL <- LANL_seq_dna[rownames(LANL_seq_dna) %in% match$LANL_ID,]

  ### Save the tmp result
  name_tmp <- paste(path_output,gsub(x=basename(name),".fasta","_closest.fasta"),sep="")
  ape::write.dna(closest_LANL,name_tmp,format="fa",colsep='', nbcol=-1)

  ### Increase the number of closest sequences if there is not enough LANL sequences
  closest_LANL <- ape::read.dna(name_tmp,format="fa")
  new_seq      <- ape::read.dna(Local_seq[i],format="fa")
  tmp <- round(nrow(closest_LANL) / (length(new_seq)+nrow(closest_LANL)) * 100,2)
  if(!is.null(nbre_close)){
   if(tmp < min_percent){
    if(verbose) cat("\t\tThere are not enough LANL sequences. \n\t\tTry with a larger nbre_close. \n")
    for(n in seq(nbre_close+10,max_closest,by=10)){
     closest  <- data.table::as.data.table(
      utils::read.csv(gsub(x=name,pattern=".fasta","_closest.txt"),
                      header=F,sep="\t",stringsAsFactors=FALSE))
     data.table::setnames(closest, paste('V',1:12,sep=''),
                          c('new_seq_ID','LANL_ID','pident','length','mismatch','gapopen',
                            'qstart','qend','sstart','send','evalue','bitscore'))
     index <- rep(1:n,length(unique(closest$new_seq_ID))) +
      rep(1:length(unique(closest$new_seq_ID))-1,each=n) * max_closest
     closest<- closest[index, ]

     match        <- base::subset(closest,select=c(new_seq_ID,LANL_ID))
     closest_LANL <- LANL_seq_dna[rownames(LANL_seq_dna) %in% match$LANL_ID,]
     name_tmp <- paste(path_output,gsub(x=basename(name),".fasta","_closest.fasta"),sep="")
     ape::write.dna(closest_LANL,name_tmp,format="fa",colsep='', nbcol=-1)

     tmp2 <- round(nrow(closest_LANL) / (length(new_seq)+nrow(closest_LANL)) * 100,2)
     if(tmp2 >= min_percent) break
    }
    if(verbose) cat("\t\tFinally keep the ",n," closest sequences.\n",sep="")
    nbre_close_tmp <- n
   }
  }

  ### Check if HXB2 is in the file
  if(!is.null(name_HXB2_file)){
   HXB2 <- ape::read.dna(name_HXB2_file,format='fa')
   name_HXB2 <- rownames(HXB2)

   if(!(name_HXB2 %in% rownames(closest_LANL))){ # Add the HXB2 sequence in the file
    if(verbose) cat("\tAdd the HXB2 sequence in the file.\n")
    name_tmp2 <- gsub(x=name_tmp,".fasta","_tmp.fasta")
    cmd <- paste("cat ",name_HXB2_file," ",name_tmp," > ",name_tmp2,'\n ',sep="")
    cmd <- paste(cmd,"rm ",name_tmp,"\n ",sep="")
    cmd <- paste(cmd,"mv ",name_tmp2," ",name_tmp,"\n ",sep="")
    system(cmd)
   }
  }

  ### Save the results
  name_match_file   <- paste(path_output,gsub(x=basename(name),".fasta","_match.txt"),sep="")
  utils::write.table(match,name_match_file,row.names=F)


  # name_result <- paste(path_output,prefixe_newseq,"_",basename(name),sep="")
  name_result <- paste(path_output,output_names[i],".fasta",sep="")
  cmd <- paste("cat ",name_tmp," ",Local_seq[i]," > ",name_result,sep="")
  system(cmd)

  # Save the LANL sequences proportion
  prop <- c(prop,tmp)
  nbre_close_res <- c(nbre_close_res,nbre_close_tmp)
 }
 return(list(prop       = prop,
             nbre_close = nbre_close_res))
}
############ ---- add_outgroup ---- ############
#' @title add_outgroup
#' @description Script for:
#' @return dqsfds
#' @param subtypes XXXX
#' @param path_output XXXX
#' @param LANL_seq XXXX
#' @param output_names XXXX#'
#' @param outgroup_size XXXX
#' @param index XXXX
#' @param verbose XXXX
#' @export
add_outgroup <- function(subtypes,
                         path_output,LANL_seq,output_names,
                         outgroup_size = 5,
                         index=NA,
                         verbose=TRUE){
 if(index == "NA") index <- NA

 if(verbose) cat("Add outgroups to each sequence files.\n")

 load(system.file(package="phyloHIV","extdata","nucl_div.rda"))

 ### Stop if only one subtype
 stopifnot(length(subtypes)>1)

 ### Stop if "subtypes" and "LANL_seq" do not have the same length
 stopifnot(length(subtypes) == length(LANL_seq))

 ### Which files ?
 if(is.na(index)){
  index <- 1:length(subtypes)
 }

 ### Define the result which will contains the name of added sequences
 res <- list()
 length(res) <- length(index)

 for(j in 1:length(index)){
  i <- index[j]
  if(i > length(subtypes)){
   stop("The index job is too large.\n")
  }

  ### Find the closest subtype
  if(length(subtypes) == 2){
   closest_subtype <- subtypes[-i]
  }else{
   tmp <- nucl_div[subtypes[-i],subtypes[i]]
   closest_subtype <- names(which.min(tmp))
  }
  index_closest <- which(subtypes == closest_subtype)

  ### Give a name for the result
  names(res)[j]  <- paste(subtypes[i],"_from_",closest_subtype,sep="")

  ### Load closest subtype sequences
  seq_to_add <- ape::read.dna(LANL_seq[index_closest],format="fa")
  if(is.list(seq_to_add))
   seq_to_add <- seq_to_add[sample(length(seq_to_add),outgroup_size)]
  if(is.matrix(seq_to_add))
   seq_to_add <- seq_to_add[sample(nrow(seq_to_add),outgroup_size),]

  ### Add these seq to the file
  name_seq_to_add <- paste(path_output,
                           "seq_to_add_to_",names(res)[j],".fasta",sep="")
  name_file <- paste(path_output,
                     output_names[i],".fasta",sep="")
  name_file_tmp <- paste(path_output,
                         output_names[i],"_tmp.fasta",sep="")

  if(verbose) cat("\tSubtype ",subtypes[i],": Add ",outgroup_size," sequences from subtype ",
                  closest_subtype," into the file ",basename(name_file),"\n",sep="")

  ape::write.dna(seq_to_add,name_seq_to_add,format="fa",colsep='', nbcol=-1)
  cmd <- paste("cat ",name_file," ",name_seq_to_add," > ",name_file_tmp,"\n ",sep="")
  cmd <- paste(cmd,'rm ',name_file,"\n ",sep="")
  cmd <- paste(cmd,'mv ',name_file_tmp," ",name_file,"\n ",sep="")
  cmd <- paste(cmd,'rm ',name_seq_to_add,sep="")
  system(cmd)

  ### Save the name of added sequences
  if(is.list(seq_to_add))
   res[[j]] <- names(seq_to_add)
  if(is.matrix(seq_to_add))
   res[[j]] <- rownames(seq_to_add)
 }

 ### Return the result
 return(res)
}


############ ---- align_seq ---- ############
#' @title align_seq
#' @description Script for:
#' @return dqsfds
#' @param names_file XXXX
#' @param aligner XXXX
#' @param thread XXXX
#' @param index XXXX
#' @param verbose XXXX
#' @export
align_seq <- function(names_file,
                      aligner,thread=5,
                      index=NA,
                      verbose=TRUE
){
 if(index == "NA") index <- NA

 if(verbose) cat("Align the sequences with ",aligner,".\n",sep="")
 path_fasta_splitter <- system.file(package="phyloHIV","ext","fasta-splitter.pl")

 ### Which files ?
 if(is.na(index)){
  index <- 1:length(names_file)
 }

 ### For each file
 for(i in index){
  if(i > length(names_file)){
   warning("The index job is too large.\n")
  }else{
   if(verbose) cat("\tFile: ",names_file[i],"\n",sep="")

   # chose the aligner
   if(aligner == "mafft"){# Split in small files if necessary

    ### Load data
    new_and_LANL <- ape::read.dna(names_file[i],format="fa")

    ### Align
    name_alignd  <- gsub(x=names_file[i],".fasta","")

    if(is.matrix(new_and_LANL)) length_new_and_LANL <- nrow(new_and_LANL) else
     length_new_and_LANL <- length(new_and_LANL)
    if( length_new_and_LANL > 500 ){
     # Build a tempory directory
     tmpdir <- paste('tmp_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),"/",sep='')
     dir.create(tmpdir)

     # Split the big file
     n_part <- length(new_and_LANL) %/% 500 + 1

     cmd <- paste(path_fasta_splitter," --n-parts ",n_part,
                  " --part-num-prefix _tmp_ ","--out-dir ",tmpdir," ",
                  names_file[i],sep="")
     system(cmd,ignore.stderr=T)

     # Align the first split
     ### do something to detect the number of 0
     cmd <- paste("ls ",tmpdir,sep="")
     names_split <- system(cmd,intern=T)
     names_split_aligned <- gsub(x=names_split,".fasta","_aligned_mafft.fasta")

     if(verbose) cat("\t\tAlignment: 1/",length(names_split),"\n",sep="")
     cmd <- paste("mafft --auto --thread ",thread," ",tmpdir,names_split[1],
                  " > ",tmpdir,names_split_aligned,
                  sep="")
     system(cmd,ignore.stderr=T)

     for(j in 2:length(names_split)){
      if(verbose) cat("\t\tAlignment: ",j,"/",length(names_split),"\n",sep="")
      cmd <- paste("mafft --auto --thread ",thread," --add ",tmpdir,names_split[j]," ",
                   tmpdir,names_split_aligned[j-1]," > ",
                   tmpdir,names_split_aligned[j],sep="")
      system(cmd,ignore.stderr=T)
      cmd <- paste("rm ",tmpdir,names_split_aligned[j-1],sep="")
      system(cmd)
     }

     # Move the result outside the tempory directory
     cmd <- paste("mv ",tmpdir,utils::tail(names_split_aligned,1)," ",name_alignd,
                  "_mafft_aligned.fasta",sep="")
     system(cmd)

     # Remove the temporary directory
     cmd <- paste("rm -r ",tmpdir,sep="")
     system(cmd)
    }else{
     cmd <- paste("mafft --auto --thread ",thread," ",names_file[i]," > ",name_alignd,"_mafft_aligned.fasta",sep="")
    }
    system(cmd,ignore.stderr=T)
   }
   # if(aligner == "CodonAlign"){
   #  ### Load data
   #  new_and_LANL <- read.dna(names_file[i],format="fa")
   #
   #  ### Align
   #  name_alignd  <- gsub(x=names_file[i],".fasta","")
   #
   #  cmd <- paste(path_programs,"CodonAlign/ca_wrapper.py -f ",names_file[i],
   #               " -od ",dirname(name_alignd),"/ -rn ",basename(name_alignd),"_CodonAlign",sep="")
   #  system(cmd,ignore.stderr=T)
   # }
  }
 }
}
############ ---- find_identical_sequences ---- ############
#' @title find_identical_sequences
#' @description Script for:
#' @return dqsfds
#' @param names_file_aligned XXXX
#' @param path_sequences_meta XXXX
#' @param dna_region XXXX
#' @param index XXXX
#' @param verbose XXXX
#' @export
find_identical_sequences <- function(names_file_aligned,path_sequences_meta,
                                     dna_region,
                                     index=NA,
                                     verbose=TRUE){
 if(index == "NA") index <- NA

 ### Which files ?
 if(is.na(index)){
  index <- 1:length(names_file_aligned)
 }

 sequences_meta <- data.table::data.table(readRDS(path_sequences_meta))
 if(verbose) cat("Find if there are identical sequences.\n")
 for(index_tmp in index){
  name <- names_file_aligned[index_tmp]
  if(verbose) cat("\tFile: ",basename(name),"\n")

  seq <- ape::read.dna(name,format="fa")

  index_dna <- grep(dna_region,rownames(seq))
  ids <- rownames(seq)[index_dna]

  tmp <- sequences_meta[sequences_meta$seqID %in% ids , c('seqID','newnum')]
  to_remove <- NULL
  for(unewnum in unique(tmp$newnum)){
   if(sum(tmp$newnum == unewnum) != 1){
    to_remove <- c(to_remove,which(tmp$newnum == unewnum)[-1])
   }
  }
  tmp <- tmp[to_remove,]
  ids_to_remove <- which(rownames(seq) %in% tmp$seqID)
  if(length(ids_to_remove) > 0)
   seq <- seq[-ids_to_remove,]

  index_dna <- NULL
  progress <- data.table::data.table( iter=floor(seq(1,nrow(seq)-1,le=11)) ,
                                      percent=seq(0,100,le=11) )
  for(i in 1:(nrow(seq)-1)){
   if(verbose & i %in% progress$iter ) cat("\t\t",progress[progress$iter==i,"percent"][[1]],"%\n")
   for(j in (i+1):nrow(seq)){
    if(sum(as.vector(seq[i,]) != as.vector(seq[j,])) == 0)
     index_dna <- rbind(index_dna,c(i,j))
   }
  }
  if(length(index_dna) != 0){
   index_dna <- unique(index_dna[,2])
   warning("\n Some duplicate sequences for this file: \n\t",name)
   if(verbose){
    cat("\t -The following sequences are removed from the file '",
        basename(name),"'.\n",sep="")
    if(length(index_dna) < 10) cat("\t\t",index_dna,"\n",sep=" ") else
     cat("\t\t",index_dna[1:10],"... \n\t\t(",length(index_dna)-10,"are hidden. )\n",sep=" ")
   }
   seq <- seq[-index_dna,]
   # write.dna(seq,gsub(x=name,".fasta","2.fasta"),format='fasta',colsep='', nbcol=-1)
   ape::write.dna(seq, name,format='fasta',colsep='', nbcol=-1)
  }

  res <- list(index=index_dna,index_id= ids_to_remove)
 }

 return(res)
}
############ ---- rm_DRM ---- ############
#' @title rm_DRM
#' @description Script for:
#' @return dqsfds
#' @param names_file_aligned XXXX
#' @param name_HXB2 XXXX
#' @param index XXXX
#' @param verbose XXXX
#' @export
rm_DRM <- function(names_file_aligned,name_HXB2,
                   index=NA,
                   verbose=TRUE){
 if(index == "NA") index <- NA

 ### Which files ?
 if(is.na(index)){
  index <- 1:length(names_file_aligned)
 }

 if(verbose) cat("Remove the Drug Resistance Mutations.\n")

 for(i in index){
  name <- names_file_aligned[i]

  if(verbose) cat("\tFile: ",basename(name),"\n",sep="")
  seq <- ape::read.dna(name,format="fa")

  tmp <- seq_rm_drugresistance2(seq=seq,
                                name_HXB2 =name_HXB2)
  sq		<- tmp$nodr.seq
  ape::write.dna(sq, gsub(x=name,".fasta","_ndrm.fasta"),format='fasta',colsep="",nbcol=-1)

 }
}
############ ---- seq_rm_drugresistance2 ---- ############
#' @title seq_rm_drugresistance2
#' @description Script for:
#' @return dqsfds
#' @param seq XXXX
#' @param outfile XXXX
#' @param name_HXB2 XXXX
#' @export
seq_rm_drugresistance2 <- function(seq, outfile=NA,name_HXB2) #="B.FR.1983.HXB2-LAI-IIIB-BRU.K03455.M" + "B.FR.1983.HXB2_copy_LAI.AF033819.M"
{
 requireNamespace("big.phylo", quietly = TRUE)
 # stopifnot( any(rownames(seq)=='HXB2') )                                      #	expect HXB2 as reference in alignment
 if(length(name_HXB2) == 1){
  if(is.list(seq)){
   index_HXB2 <- grep(names(seq),pattern=name_HXB2)[1] # PMG 12/04/19
  }else{
   index_HXB2 <- grep(rownames(seq),pattern=name_HXB2)[1] # PMG 03/09/18    # enlever le [1] normalement
  }
 }
 if(length(name_HXB2) == 2){
  if(is.list(seq)){
   index_HXB2 <- c(grep(names(seq),pattern=name_HXB2[1]),
                   grep(names(seq),pattern=name_HXB2[2]))[1] # PMG 12/04/19
  }else{
   index_HXB2 <- c(grep(rownames(seq),pattern=name_HXB2[1]),
                   grep(rownames(seq),pattern=name_HXB2[2]))[1] # PMG 03/09/18    # enlever le [1] normalement
  }
 }
 stopifnot( length(index_HXB2) != 0 ) # PMG 03/09/18                            #	expect HXB2 as reference in alignment

 load(system.file(package="big.phylo", 'AC_drugresistance_201508.rda'))
 load(system.file(package="big.phylo", 'refseq_hiv1_hxb2.rda'))

 index_HXB2.K03455 <- which(names(hxb2) == "HXB2.K03455")
 hxb2			  <- paste(hxb2[, index_HXB2.K03455], collapse='')
 if(is.list(seq)){
  seq.hxb2	<- paste(as.character(seq[index_HXB2]),collapse='')
 }else{
  seq.hxb2	<- paste(as.character(seq[index_HXB2,]),collapse='')
 }

 #	check that HXB2 in alignment is HXB2
 seq.hxb2.ng	<- gsub('-+','',seq.hxb2)
 seq.hxb2.ng <- gsub('n' ,'',seq.hxb2.ng) # PMG 04/09/18
 stopifnot(nchar(seq.hxb2.ng)<=nchar(hxb2))					                                #	expect HXB2 in alignment is part of true HXB2
 seq.hxb2.st		<- as.integer(regexpr(substr(seq.hxb2.ng,1,20),hxb2))
 stopifnot(seq.hxb2.st>0)						                                                	#	expect HXB2 in alignment is part of true HXB2
 stopifnot(nchar(seq.hxb2.ng)+seq.hxb2.st-1L<=nchar(hxb2))	                     # expect HXB2 in alignment is part of true HXB2
 stopifnot( seq.hxb2.ng==substr(hxb2, seq.hxb2.st, nchar(seq.hxb2.ng)+seq.hxb2.st-1L) )

 #	find coordinates of real HXB2 in alignment
 seq.hxb2.pos	<- c()
 tmp				<- gregexpr('-+',seq.hxb2)[[1]]
 if(tmp[1]>-1)
 {
  seq.hxb2.pos <- data.table::data.table(GP_ST = as.integer(tmp),
                             GP_LEN = attr(tmp,'match.length'))
  seq.hxb2.pos <- seq.hxb2.pos[, list(GP_POS = seq.int(GP_ST, length.out = GP_LEN)),
                               by='GP_ST'][, GP_POS]
 }

 seq.hxb2.pos 	<- data.table::data.table(HXB2INSEQ_POS = setdiff(seq_len(nchar(seq.hxb2)),
                                                     seq.hxb2.pos))
 # seq.hxb2.pos[, HXB2_POS := seq_len(nrow(seq.hxb2.pos)) ]
 seq.hxb2.pos$HXB2_POS <- seq_len(nrow(seq.hxb2.pos))

 data.table::set(seq.hxb2.pos, NULL, 'HXB2_POS',
                 seq.hxb2.pos[, "HXB2_POS"]+seq.hxb2.st-1L )
 stopifnot( utils::tail(seq.hxb2.pos[,"HXB2_POS" ],1)<=nchar(hxb2) )
 # stopifnot( seq.hxb2.pos[, tail(HXB2_POS,1)]<=nchar(hxb2) )


 #	merge with dr
 data.table::setnames(dr, 'HXB2.pos', 'HXB2_POS')
 dr  <- merge(dr, seq.hxb2.pos, by='HXB2_POS')
 data.table::setnames(dr, 'HXB2INSEQ_POS', 'Alignment.nuc.pos')
 utils::capture.output(tmp <- big.phylo:::seq.rm.drugresistance.internal(as.character(seq), dr, verbose=F, rtn.DNAbin=1 )) ### XXX
 # capture.output is used above in order to avoid any print # PMG 10-09-2018
 dr.info		<- tmp$nodr.info
 nodr.seq	<- tmp$nodr.seq
 #
 #	save
 #
 if(!is.na(outfile))
  save(seq, nodr.seq, dr.info, file=outfile)
 list(nodr.info=dr.info, nodr.seq=nodr.seq)
}
