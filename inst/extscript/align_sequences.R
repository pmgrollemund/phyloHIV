########################################################################### ----
##################### Load Rscript options #####################################
########################################################################### ----
require(rjson)
require(phyloHIV)
require(big.phylo)

# args <- commandArgs(trailingOnly=TRUE)
# args <- c("./extjson/options.json")
args <- c("./json/options.json")
if(length(args)>0){
 json_file <- args[1]
 json_data <- fromJSON(file=json_file)

 for(i in 1:length(json_data)){
  assign(names(json_data)[i], json_data[[i]])
 }
}

########################################################################### ----
##################### Initialization ###########################################
########################################################################### ----
package_name <- "phyloHIV"

### Load internal options ----
json_data <- fromJSON(file=system.file(package=package_name,"extjson","internal_options.json"))
for(i in 1:length(json_data)){
 assign(names(json_data)[i], json_data[[i]])
}

### Redefine some names and paths ----
LANL_seq  <- paste(path_LANL_seq,LANL_seq_names,sep="")
LANL_db <- paste(path_LANL_seq,
                 gsub(x=LANL_seq_names,
                      pattern=".fasta",
                      replacement=".db"),sep="")
path_sequences_meta <- paste(path_RDS,RDS_names[1],sep="")

Local_seq <- paste(path_local_seq,local_seq_names,sep="")

name_HXB2_file <- system.file(package=package_name,"extdata","LANL","HXB2.fasta")

########################################################################### ----
##################### Get the closest sequences ################################
########################################################################### ----

### Make a blast database ----
makeblastdb(LANL_seq,LANL_db,
            job_index,verbose)
# Run the following command for each subtype: (mixing between R and Linux command)
#     makeblastdb -in inputs -dbtype nucl -title
#       gsub(".fasta","",inputs)
#     -out outputs
# Where:
#     -in: the sequence file to make the database
#     -dbtype: the type of the sequences (nucleotide)
#     -title : the title of the job
#     -out : the name of the output file

### Run blastn ----
blastn(Local_seq,LANL_db,max_closest,
       job_index,verbose)
# Run the following command for each subtype: (mixing between R and Linux command)
#     blastn -query name_Subtype_new_sequences -db names_db -out
#       gsub("_stripped.db","_closest.txt",names_db)
#     -max_target_seqs nbre_close -outfmt 6
# Where:
#     -query: the sequence files ("studied sequences") for which we are looking for the closest in names_db
#     -db: the blast database
#     -out: the name of the output
#     -max_target_seqs: the number of closest sequences to find in the database for each studied sequence

### Get closest sequences ----
closest_result <- get_closest_seq(LANL_seq,Local_seq,
                                  path_aligned,seq_names,
                                  max_closest, # maximal number of closest sequences
                                  nbre_close,   # number of closest sequences
                                  min_percent, # minimal percent of closest sequences
                                  name_HXB2_file,
                                  job_index,verbose)
if(is.na(job_index) || job_index == "NA"){
 save(closest_result,file=paste(path_aligned,"Closest_result.rda",sep=""))
}else{
 save(closest_result,file=paste(path_aligned,"Closest_result_",job_index,".rda",sep=""))
}

### Add an outgroup sequences from the closest subtype ----
outgroup <- add_outgroup(subtypes,
                         path_aligned,LANL_seq,seq_names,
                         outgroup_size,
                         job_index,verbose)
if(is.na(job_index) || job_index == "NA"){
 save(outgroup,file=paste(path_aligned,"outgroup.rda",sep=""))
}else{
 save(outgroup,file=paste(path_aligned,"outgroup_",job_index,".rda",sep=""))
}

### Check that HXB2 is in each file  ----
names_file <- paste(path_aligned,seq_names,".fasta",sep="")
check_HXB2(names_file,name_HXB2_file,
           job_index,verbose)

########################################################################### ----
##################### Align the sequences ######################################
########################################################################### ----
### Align ----
align_seq(names_file,aligner,thread,
          job_index,verbose)
# Run the following command for each subtype: (mixing between R and Linux command)
# mafft --auto --thread 5
#      --add tmpdir/names_split[j]
#      tmpdir/names_split_aligned[j-1] > tmpdir/names_split_aligned[j]
# Where:
#     --auto : to get the defaut options in terms of speed and accuracy
#     --thread: number of threads (parallele running)
#     -add : to merge the output to a given file

### Remove some gaps column ----
names_file_aligned <- gsub(x=names_file,".fasta",paste("_",aligner,"_aligned.fasta",sep=""))

for(name in names_file_aligned){
 seq <- read.dna(name,format="fa")
 seq	<- seq.strip.gap(seq, strip.pc=1-nbre_col_rm_strip/nrow(seq))
 write.dna(seq,name,format='fasta',colsep='', nbcol=-1)
}
### Remove more columns ----
# for(name in names_file_aligned[4]){
#  seq <- read.dna(name,format="fa")
#  seq  <- seq.strip.gap(seq, strip.pc=1-100/nrow(seq))
#  write.dna(seq,name,format='fasta',colsep='', nbcol=-1)
# }

### Remove the identical sequences ----
indexes <- find_identical_sequences(names_file_aligned,path_sequences_meta,
                                    dna_region,
                                    job_index,verbose)
if(is.na(job_index)){
 save(indexes,file=paste(path_aligned,"indexes_identical.rda",sep=""))
}else{
 save(indexes,file=paste(path_aligned,"indexes_identical_",job_index,".rda",sep=""))
}

# ### Remove Drug Resistance Mutations ----
rm_DRM(names_file_aligned,name_HXB2,
       job_index,verbose)
