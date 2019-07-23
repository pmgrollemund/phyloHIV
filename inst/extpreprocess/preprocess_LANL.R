########################################################################### ----
##################### Load Rscript options #####################################
########################################################################### ----
require(rjson)
require(phyloHIV)
require(big.phylo)

# args <- commandArgs(trailingOnly=TRUE)
# args <- c("./extjson/options.json")
# args <- system.file(package="phyloHIV","extpreprocess","options_preprocess.json")
args <- "/home/paul/Documents/Postdoc/HIV_transmission/Code/phyloHIV/inst/extpreprocess/options_preprocess.json"
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
### Redefine some names and paths ----
path_LANL_seq <- paste(normalizePath(path_LANL_seq),"/",sep="")
LANL_seq <- paste(path_LANL_seq,LANL_seq_names,sep="")

########################################################################### ----
##################### Preprocess LANL sequences ################################
########################################################################### ----
### Remove the duplicate sequences ----
for(name_file in LANL_seq){
 clean_LANL(name_file,name_length)
}

### Remove gaps from the database sequences ----
rm_gaps_LANL(LANL_seq)
