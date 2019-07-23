############ ---- check_HXB2 ---- ############
#' @title check_HXB2
#' @description Script for:
#' @return dqsfds
#' @param name_file XXXX
#' @param name_HXB2_file XXXX
#' @param index XXXX
#' @param verbose XXXX
#' @export
check_HXB2 <- function(name_file,name_HXB2_file,
                       index=NA,
                       verbose=TRUE){
 if(index == "NA") index <- NA

 if(verbose) cat("Check that HXB2 is in each file.\n")

 ### Which files ?
 if(!is.na(index)){
  name_file <- name_file[index]
 }

 for(name in name_file){
  ### Load data
  seq  <- ape::read.dna(name,format="fa")
  HXB2 <- ape::read.dna(name_HXB2_file,format="fa")

  ### Check and add HXB2 if necessary
  name_HXB2 <- rownames(HXB2)
  if(name_HXB2 %in% names(seq)){
   if(verbose) cat("\t",basename(name)," already contains HXB2.\n",sep="")
  }else{
   if(verbose) cat("\tAdd HXB2 into the file ",basename(name),"\n",sep="")
   name_tmp <- gsub(x=name,".fasta","_tmp.fasta")

   cmd <- paste("cat ",name_HXB2_file," ",name," > ",name_tmp,"\n ",sep="")
   cmd <- paste(cmd,"rm ",name,"\n ",sep="")
   cmd <- paste(cmd,"mv ",name_tmp," ",name,"\n",sep="")
   system(cmd)
  }
 }
}
############ ---- extract_HXB2 ---- ############
#' @title extract_HXB2
#' @description Extract the HXB2 sequence from a LANL dataset.
#' @return dsfdsf
#' @param path_LANL_sequences character string giving the path to the LANL datasets.
#' @param prefixe_LANL character string giving the first characters of the LANL datasets.
#' @param name_HXB2 character string giving the name of the HXB2 sequence.
#' @param name_HXB2_file character string giving the path name of the output.
#' @param verbose Write stuff. [default = TRUE]
#' @export
extract_HXB2 <- function(path_LANL_sequences,prefixe_LANL,
                         name_HXB2,name_HXB2_file,
                         verbose=TRUE){
 if(length(name_HXB2) > 1) name_HXB2 <- name_HXB2[1]

 if(verbose) cat("Search the HXB2 sequence.\n")
 ### Find the LANL datasets.
 cmd <- paste("ls ",path_LANL_sequences,prefixe_LANL,"*.fasta",sep="")
 names <- system(cmd,intern=T)

 ### For each LANL dataset: search the HXB2 sequence.
 for(name in names){
  seq <- ape::read.dna(name,format='fa')
  if(name_HXB2 %in% rownames(seq)){
   ### If this LANL dataset contains the HXB2 sequence.
   index <- which(grepl(x=rownames(seq),name_HXB2))
   HXB2  <- seq[index,]
   if(verbose) cat("\tThe HXB2 sequence is found in the file ",name,"\n",sep="")

   ### Do not search the sequence in the other LANL datasets.
   break
  }
 }

 ### Save the HXB2 sequence in a fasta file.
 ape::write.dna(HXB2,name_HXB2_file,format='fasta',colsep='', nbcol=-1)
}
############ ---- write_qsub ---- ############
#' @title write_qsub
#' @description Extract the HXB2 sequence from a LANL dataset.
#' @return dsfdsf
#' @param Rscript XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param options_file XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param path_output XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param job_title XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param hpc.q XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param hpc.select XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param hpc.walltime XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param hpc.mem XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param hpc.nproc XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param hpc.load XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param hpc.array XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param log XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param multiple_log XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @export
write_qsub <- function(Rscript,options_file=NA,
                       path_output="./",job_title=NA,
                       hpc.q="pqeelab",hpc.select=1,hpc.walltime=24,hpc.mem="1gb",
                       hpc.nproc=1,hpc.load=big.phylo:::HPC.CX1.IMPERIAL.LOAD,
                       hpc.array=NA,
                       log=TRUE,multiple_log=TRUE
                       ){
 # Reference time
 time_tmp <- format(Sys.time(),"%y-%m-%d-%H-%M-%S")

 # The log file name
 log_name <- paste(path_output,"log_",basename(Rscript),"_",time_tmp,sep="")

 # Call a big.phylo function that writes the PBS header
 qsub_file  <- big.phylo::cmd.hpcwrapper.cx1.ic.ac.uk(
  hpc.select = hpc.select,hpc.walltime=hpc.walltime,hpc.q=hpc.q,hpc.mem=hpc.mem,
  hpc.nproc=hpc.nproc,hpc.load=NULL,hpc.array=hpc.array)

 if(!is.na(job_title)){
  qsub_file <- paste(qsub_file,"#PBS -N '",job_title,"'\n",sep="")
 }
 if(log || multiple_log){
  qsub_file <- paste(qsub_file,"#PBS -o '",log_name,".txt'\n",sep="")
 }
 if(!is.na(hpc.load)){
  qsub_file <- paste(qsub_file,"\n", hpc.load,"\n",sep = "")
 }

 # Add the Rscript call
 qsub_file <- paste(qsub_file,"\nRscript ",Rscript,sep="")
 if(!is.na(options_file)){
  qsub_file <- paste(qsub_file," ",options_file,sep="")
 }
 if(!is.na(hpc.array)){
  qsub_file <- paste(qsub_file," $PBS_ARRAY_INDEX",sep="")
  if(multiple_log){
   qsub_file <- paste(qsub_file," >> ",log_name,"_${PBS_ARRAY_INDEX}.txt",sep="")
  }
 }

 # Write the qsub file
 outfile_name <- paste("script_",basename(Rscript),"_",time_tmp,".qsub",sep="")
 write(qsub_file,file=paste(path_output,outfile_name,sep=""))
}
############ ---- write_option_file ---- ############
#' @title write_option_file
#' @description Extract the HXB2 sequence from a LANL dataset.
#' @return dsfdsf
#' @param output_file XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param path_RDS XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param RDS_names XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param path_LANL_seq XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param LANL_seq_names XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param path_seq XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param seq_names XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param path_file_aligned XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param names_file_aligned XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param path_tree XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param tree_names XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param path_colored_tree XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param path_clades XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param local_group XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param focus_group XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param bs.n XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param subtypes XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param job_index XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @param verbose XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
#' @export
write_option_file <- function(output_file = "./options.json",
                              path_RDS="./extdata/FakeData/",
                              RDS_names=c("sequences_meta_modified.rds",
                                          "person_modified.rds",
                                          "Country_db.rds"),
                              path_LANL_seq="./extdata/LANL/",
                              LANL_seq_names=c("LANL_Subtype_01_AE_stripped.fasta",
                                               "LANL_Subtype_A1_stripped.fasta"),
                              path_seq="./extdata/Subtype/",
                              seq_names=c("newseq_Subtype_01_AE.fasta",
                                          "newseq_Subtype_A1.fasta"),
                              path_file_aligned="./extdata/Closest/",
                              names_file_aligned=c("newseq_LANL_Subtype_01_AE",
                                                   "newseq_LANL_Subtype_A1"),
                              path_tree="./extdata/Tree/",
                              tree_names=c("newseq_LANL_Subtype_01_AE_mafft_aligned_ndrm",
                                          "newseq_LANL_Subtype_A1_mafft_aligned_ndrm"),
                              path_colored_tree="./extdata/Colored_tree/",
                              path_clades="./extdata/Clades/",
                              local_group="WA",
                              focus_group="KC",
                              bs.n=2,
                              subtypes=c("01_AE","A1"),
                              job_index=NA,
                              verbose="true"){
 txt <- '{\n'

 txt <- paste(txt,'"path_RDS": "',path_RDS,'",\n',sep='')
 txt <- paste(txt,'"RDS_names": ["',
              paste(RDS_names,collapse='","')
              ,'"], \n',sep='')

 txt <- paste(txt,'"path_LANL_seq": "',path_LANL_seq,'",\n',sep='')
 txt <- paste(txt,'"LANL_seq_names": ["',
              paste(LANL_seq_names,collapse='","')
              ,'"], \n',sep='')

 txt <- paste(txt,'"path_seq": "',path_seq,'",\n',sep='')
 txt <- paste(txt,'"seq_names": ["',
              paste(seq_names,collapse='","')
              ,'"], \n',sep='')

 txt <- paste(txt,'"path_file_aligned": "',path_file_aligned,'",\n',sep='')
 txt <- paste(txt,'"names_file_aligned": ["',
              paste(names_file_aligned,collapse='","')
              ,'"], \n',sep='')

 txt <- paste(txt,'"path_tree": "',path_tree,'",\n',sep='')
 txt <- paste(txt,'"tree_names": ["',
              paste(tree_names,collapse='","')
              ,'"], \n',sep='')

 txt <- paste(txt,'"path_colored_tree": "',path_colored_tree,'",\n',sep='')
 txt <- paste(txt,'"path_clades": "',path_clades,'",\n',sep='')

 txt <- paste(txt,'"local_group": "',local_group,'",\n',sep='')
 txt <- paste(txt,'"focus_group": "',focus_group,'",\n',sep='')
 txt <- paste(txt,'"bs.n": ',bs.n,',\n',sep='')
 txt <- paste(txt,'"subtypes": ["',
              paste(subtypes,collapse='","')
              ,'"], \n',sep='')
 txt <- paste(txt,'"job_index": "',job_index,'",\n',sep='')
 txt <- paste(txt,'"verbose": ',verbose,',\n',sep='')

 txt <- paste(txt,'}',sep='')

 # Write the option file
 cat(txt,file=output_file)
}


############ ---- clean_LANL ---- ############
#' @title Clean the LANL database.
#' @description Remove sequences which are duplicates (with regards to accesion
#'   number or name) and remove sequences for which the name does not contain
#'   the required information.
#' @param name_file a character vector, the LANL file names by subtype.
#' @param info_number an integer, the number of information in the sequence names.
#' @param verbose Write stuff.
#' @export
clean_LANL <- function(name_file,info_number,
                       verbose = TRUE){
 ### Save data
 if(verbose) cat(paste(basename(name_file),": \n",sep=""))
 if(verbose) cat("\t Save original version of the LANL data.\n")
 cmd <- paste("cp ",name_file," ",
              paste(gsub(pattern=".fas|.fasta",replacement='',name_file),
                    "_save.fasta",sep=""))
 system(cmd)

 ### Load data
 if(verbose) cat("\t Clean the LANL dataset.\n")
 LANL_sequences <- ape::read.dna(name_file,format="fa")

 ### Remove sequences whose the name is too long
 index <- which( (unlist(lapply(strsplit(rownames(LANL_sequences),split="[.]"),length)))
                 != info_number )
 if(length(index) != 0 ){
  LANL_sequences <- LANL_sequences[-index,]
 }

 ### Remove sequences which are about the sames individuals (accession number)
 accession_number <- unlist( lapply( strsplit(rownames(LANL_sequences),split="[.]") , function(vec) vec[5]) )

 # If accesion_number is not unique: find the duplicates
 if(!(is.unique(accession_number))){
  index <- NULL
  for(i in 1:(length(accession_number)-1)){
   tmp <- which(accession_number == accession_number[i])
   if( length(tmp) != 1 ){
    tmp <- tmp[tmp>i]
    index <- c(index,tmp)
   }
  }
  index <- unique(index)

  # Remove the duplicates
  LANL_sequences <- LANL_sequences[-index,]
 }

 ## Remove sequences which are about the sames individuals (name)
 name <- unlist( lapply( strsplit(rownames(LANL_sequences),split="[.]") , function(vec) vec[4]) )

 # If name is not unique: find the duplicates
 if(!(is.unique(name))){
  index <- NULL
  for(i in 1:(length(name)-1)){
   tmp <- which(name == name[i])
   if( length(tmp) != 1 ){
    tmp <- tmp[tmp>i]
    index <- c(index,tmp)
   }
  }
  index <- unique(index)

  # Remove the duplicates
  LANL_sequences <- LANL_sequences[-index,]
 }

 ### Save result
 ape::write.dna(LANL_sequences,name_file,format="fasta",colsep='', nbcol=-1)

}

############ ---- is.unique ---- ############
#' @title is.unique
#' @description Script for:
#' \describe{
#'   \item{\code{XXX}}{XXX}
#'   }
#' @param v XXXX
#' @export
is.unique <- function(v) length(unique(v)) == length(v)


############ ---- rm_gaps_LANL ---- ############
#' @title rm_gaps_LANL
#' @description Script for:
#' \describe{
#'   \item{\code{XXX}}{XXX}
#'   }
#' @param LANL_seq XXXX
#' @param code.from XXXX
#' @param code.to XXXX
#' @param verbose XXXX
#' @export
rm_gaps_LANL <- function(LANL_seq,code.from="-",code.to="n",
                         verbose=TRUE){
 if(verbose) cat("Remove gaps from the database sequences.\n")

 ### For each subtype
 for(i in 1:length(LANL_seq)){
  name_file <- LANL_seq[i]
  names_result <- name_file
  if(verbose) cat(paste("\tFile: ",name_file,"\n",sep=""))
  ### Load data
  tmp <- ape::read.dna(name_file,format="fa")
  tmp	<- big.phylo::seq.strip.gap(tmp, strip.pc=1-10/nrow(tmp))
  tmp	<- big.phylo::seq.replace(tmp, code.from=code.from, code.to=code.to, verbose=0)
  ape::write.dna(tmp,file=names_result,format='fasta',colsep='', nbcol=-1)
 }
}
