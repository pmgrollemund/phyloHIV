########################################################################### ----
##################### Load Rscript options #####################################
########################################################################### ----
require(rjson)
require(phyloHIV)
require(big.phylo)

# args <- c("inst/tree_options.json",1)

args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
 json_file <- args[1]
 json_data <- fromJSON(file=json_file)

 for(i in 1:length(json_data)){
  assign(names(json_data)[i], json_data[[i]])
 }
}

if(length(args) > 1){
 index <- as.numeric(args[2])
}else{
 index <- NA
}

########################################################################### ----
##################### Initialization ###########################################
########################################################################### ----
package_name <- "phyloHIV"

### Load internal options ----
json_data <- fromJSON(file=system.file(package=package_name,"internal_options.json"))
for(i in 1:length(json_data)){
 assign(names(json_data)[i], json_data[[i]])
}

### Define names and paths ----
names_file_aligned_ndrm <- paste(path_file_aligned,names_file_aligned,sep="")
path_outgroup <- paste(path_file_aligned,"outgroup.rda",sep="")

### Which file to treat ----
all_file_bs <- expand.grid(names_file_aligned_ndrm,1:bs.n)
names(all_file_bs) <- c("name","bs.id")

if(!is.na(index)){
 names_file_aligned_ndrm <- as.character(all_file_bs[index,1])
 bs.id <- as.numeric(all_file_bs[index,2])
}else{
 bs.id <- NA
}

########################################################################### ----
##################### Build the trees ##########################################
########################################################################### ----
### Building trees ----
building_tree(names_file_aligned_ndrm,path_output,
              tree_builder,
              bs.n,bs.id,
              verbose)

### Postprocess trees (rooting and remove outgroup sequences) ----
tree_names <-
 paste(path_output,"fasttree_",gsub(x=basename(names_file_aligned_ndrm),".fasta","/"),
       gsub(x=basename(names_file_aligned_ndrm),".fasta",".newick"),sep='')

postprocess_tree(tree_names,path_outgroup,name_HXB2,
                 bs.n,bs.id,
                 verbose)
