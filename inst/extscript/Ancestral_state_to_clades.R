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

### Which file to treat ----
all_file_bs <- expand.grid(tree_names,1:bs.n)
names(all_file_bs) <- c("tree_names","bs.id")
all_file_bs$tree_names <- paste(normalizePath(path_tree),"/",tree_builder,"_",
                                all_file_bs$tree_names,"/",all_file_bs$tree_names,
                                "_rerooted.newick",sep="")
all_file_bs$renamed_tree_names <- gsub(x=all_file_bs$tree_names,
                                       pattern=".newick",
                                       replacement="_renamed.newick")

if(!is.na(job_index) && job_index != "NA"){
 tree_names <- as.character(all_file_bs$tree_names[job_index])
 bs.id <- as.numeric(all_file_bs$bs.id[job_index])
 renamed_tree_names <- as.character(all_file_bs$renamed_tree_names[job_index])
}else{
 tree_names <- unique(as.character(all_file_bs$tree_names))
 bs.id <- NA
 renamed_tree_names <- unique(as.character(all_file_bs$renamed_tree_names))
}

### Load and define objects ----
path_RDS <- paste(normalizePath(path_RDS),"/",sep="")
sequences_meta <- data.table(readRDS(file.path(path_RDS,RDS_names[1])))
person         <- data.table(readRDS(file.path(path_RDS,RDS_names[2])))
Country_db     <- readRDS(file.path(path_RDS,RDS_names[3]))

########################################################################### ----
##################### Reconstruct the ancestral states #########################
########################################################################### ----
#- Begin relabelling -#
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
#- End relabelling -#

#- Begin phyloscanner -#
### Phyloscanner ----
wrap_phyloscanner_analyse_tree(renamed_tree_names,path_colored_tree,
                               opt_regex=opt_regex,
                               bs.n = bs.n,bs.id = bs.id,
                               verbose =verbose)
#- End phyloscanner -#
path_colored_tree <- paste(normalizePath(path_colored_tree),"/",sep="")

#- Begin postprocess_phyloscanner -#
### Phyloscanner postprocess ----
if(!is.null(focus_subgroup)) focus_group <- paste(focus_group,focus_subgroup,sep="")
phyloscanner_postprocess(renamed_tree_names,path_colored_tree,
                         focus_group=focus_group,separator=separator,
                         bs.n = bs.n,bs.id = bs.id,
                         verbose =verbose)
#- End postprocess_phyloscanner -#

########################################################################### ----
##################### Extract the clades ######################################
########################################################################### ----
#- Begin extract_clades -#
### Clades extracting for each bootstrap index ----
clades_extracting(renamed_tree_names,
                  path_colored_tree,path_clades,
                  focus_group,
                  bs.n=bs.n,bs.id=bs.id,
                  verbose=verbose)
#- End extract_clades -#
