########################################################################### ----
##################### Load Rscript options #####################################
########################################################################### ----
require(rjson)
require(phyloHIV)
require(big.phylo)

args <- commandArgs(trailingOnly=TRUE)
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

### Define names and paths ----
names_file_aligned_ndrm <- paste(normalizePath(path_aligned),"/",tree_names,".fasta",sep="")
path_outgroup <- paste(normalizePath(path_aligned),"/","outgroup.rda",sep="")

### Which file to treat ----
all_file_bs <- expand.grid(names_file_aligned_ndrm,1:bs.n)
names(all_file_bs) <- c("name","bs.id")

if(!is.na(job_index) && job_index != "NA"){
 names_file_aligned_ndrm <- as.character(all_file_bs[job_index,1])
 bs.id <- as.numeric(all_file_bs[job_index,2])
}else{
 bs.id <- NA
}

########################################################################### ----
##################### Build the trees ##########################################
########################################################################### ----
#- Begin build_tree -#
### Building trees ----
building_tree(names_file_aligned_ndrm,path_tree,
              tree_builder,
              bs.n,bs.id,
              verbose)
path_tree <- paste(normalizePath(path_tree),"/",sep="")
#- End build_tree -#

#- Begin postprocess_tree -#
### Postprocess trees (rooting and remove outgroup sequences) ----
tree_names <-
 paste(path_tree,"fasttree_",gsub(x=basename(names_file_aligned_ndrm),".fasta","/"),
       gsub(x=basename(names_file_aligned_ndrm),".fasta",".newick"),sep='')

postprocess_tree(tree_names,path_outgroup,name_HXB2,
                 bs.n,bs.id,
                 verbose)
#- End postprocess_tree -#
