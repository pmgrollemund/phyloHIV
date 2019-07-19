########################################################################### ----
##################### Load Rscript options #####################################
########################################################################### ----
require(rjson)
require(phyloHIV)
require(big.phylo)

# args <- c("inst/clades_options.json")

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

### Which file to treat ----
all_file_bs <- expand.grid(tree_names,1:bs.n)
names(all_file_bs) <- c("tree_names","bs.id")
all_file_bs$tree_names <- paste(path_tree,tree_builder,"_",
                                all_file_bs$tree_names,"/",all_file_bs$tree_names,
                                "_rerooted.newick",sep="")
all_file_bs$renamed_tree_names <- gsub(x=all_file_bs$tree_names,
                                       pattern=".newick",
                                       replacement="_renamed.newick")

if(!is.na(index)){
 tree_names <- as.character(all_file_bs$tree_names[index])
 bs.id <- as.numeric(all_file_bs$bs.id[index])
 renamed_tree_names <- as.character(all_file_bs$renamed_tree_names[index])
}else{
 tree_names <- unique(as.character(all_file_bs$tree_names))
 bs.id <- NA
 renamed_tree_names <- unique(as.character(all_file_bs$renamed_tree_names))
}

### Load and define objects ----
sequences_meta <- data.table(readRDS(file.path(path_RDS,RDS_names[1])))
person         <- data.table(readRDS(file.path(path_RDS,RDS_names[2])))
Country_db     <- readRDS(file.path(path_RDS,RDS_names[3]))

########################################################################### ----
##################### Reconstruct the ancestral states #########################
########################################################################### ----
### Relabelling ----
tips_relabelling(tree_names,
                 person,sequences_meta,Country_db,
                 vars_for_label,
                 local_group=local_group,
                 focus_group=focus_group,
                 focus_subgroup=focus_subgroup,
                 index_subgroup=index_subgroup,
                 name_subgroup=name_subgroup,
                 negative=negative,
                 local_indicator=local_indicator,
                 focus_indicator=focus_indicator,
                 separator=separator,
                 blacklist_paths=blacklist_paths,
                 blacklist_names=blacklist_names,
                 Unknown = Unknown,
                 local_tip_split=local_tip_split,
                 external_tip_split=external_tip_split,
                 local_index_location=local_index_location,
                 external_index_location=external_index_location,
                 bs.n=bs.n,bs.id=bs.id,
                 verbose =verbose)

### Phyloscanner ----
wrap_phyloscanner_analyse_tree(renamed_tree_names,path_colored_tree,
                               opt_regex=opt_regex,
                               bs.n = bs.n,bs.id = bs.id,
                               verbose =verbose)

### Phyloscanner postprocess ----
if(!is.null(focus_subgroup)) focus_group <- paste(focus_group,focus_subgroup,sep="")
phyloscanner_postprocess(renamed_tree_names,path_colored_tree,
                         focus_group=focus_group,
                         bs.n = bs.n,bs.id = bs.id,
                         verbose =verbose)

########################################################################### ----
##################### Extract the clades ######################################
########################################################################### ----
### Clades extracting for each bootstrap index ----
clades_extracting( unique(as.character(all_file_bs$renamed_tree_names)),
                  path_colored_tree,path_clades,
                  focus_group,
                  bs.n=bs.n,bs.id=bs.id,
                  verbose=verbose)
