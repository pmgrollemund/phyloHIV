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


# Rscript <- system.file(package="phyloHIV","extscript","align_sequences.R")
# options_file <- system.file(package="phyloHIV","extjson","align_options.json")
#
# write_qsub(Rscript, options_file,
#            path_output="inst/" ,
#            hpc.array=3)
