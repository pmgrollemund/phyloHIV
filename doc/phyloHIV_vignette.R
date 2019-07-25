## ---- include = FALSE-----------------------------------------------------------------------------
options(width = 100)
knitr::opts_chunk$set(echo = TRUE, dev = c('png','pdf'), fig.width= 6, 
                      fig.height=4,
                      collapse = TRUE,
                      comment = "#>"
)

## ---- echo=FALSE,eval=TRUE------------------------------------------------------------------------
Rscript <- readLines(file.path("..","inst","extscript","align_sequences.R"))

begin <- which(Rscript == "#- Begin blast -#") + 1 
end <- which(Rscript == "#- End blast -#") - 1
Blast_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin get_closest_seq -#") + 1 
end <- which(Rscript == "#- End get_closest_seq -#") - 1
get_closest_seq_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin align -#") + 1 
end <- which(Rscript == "#- End align -#") - 1
align_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin postprocess_align_gaps -#") + 1 
end <- which(Rscript == "#- End postprocess_align_gaps -#") - 1
postprocess_align_gaps_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin postprocess_align_identical -#") + 1 
end <- which(Rscript == "#- End postprocess_align_identical -#") - 1
postprocess_align_identical_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin postprocess_align_DRM -#") + 1 
end <- which(Rscript == "#- End postprocess_align_DRM -#") - 1
postprocess_align_DRM_Rcsript <- Rscript[begin:end]

Rscript <- readLines(file.path("..","inst","extscript","build_trees.R"))

begin <- which(Rscript == "#- Begin build_tree -#") + 1 
end <- which(Rscript == "#- End build_tree -#") - 1
build_tree_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin postprocess_tree -#") + 1 
end <- which(Rscript == "#- End postprocess_tree -#") - 1
postprocess_tree_Rcsript <- Rscript[begin:end]

Rscript <- readLines(file.path("..","inst","extscript","Ancestral_state_to_clades.R"))

begin <- which(Rscript == "#- Begin relabelling -#") + 1 
end <- which(Rscript == "#- End relabelling -#") - 1
relabelling_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin phyloscanner -#") + 1 
end <- which(Rscript == "#- End phyloscanner -#") - 1
phyloscanner_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin postprocess_phyloscanner -#") + 1 
end <- which(Rscript == "#- End postprocess_phyloscanner -#") - 1
postprocess_phyloscanner_Rcsript <- Rscript[begin:end]

begin <- which(Rscript == "#- Begin extract_clades -#") + 1 
end <- which(Rscript == "#- End extract_clades -#") - 1
extract_clades_Rcsript <- Rscript[begin:end]

## ---- echo=TRUE,eval=FALSE,code=readLines("../inst/extjson/options.json",warn=FALSE)--------------
#  {
#   "job_index": "NA",
#   "path_RDS": "./extdata/FakeData/",
#   "RDS_names": [
#                "sequences_meta_modified.rds",
#                "person_modified.rds",
#                "Country_db.rds"
#                ],
#  
#  
#   "path_LANL_seq": "./extdata/LANL/",
#   "LANL_seq_names": [
#    "LANL_Subtype_01_AE.fasta",
#    "LANL_Subtype_A1.fasta"
#    ],
#   "path_local_seq": "./extdata/Local/",
#   "local_seq_names": [
#    "newseq_Subtype_01_AE.fasta",
#    "newseq_Subtype_A1.fasta"
#    ],
#  
#  
#   "path_aligned": "./extdata/Align/",
#   "seq_names": [
#    "newseq_LANL_Subtype_01_AE",
#    "newseq_LANL_Subtype_A1"
#    ],
#  
#  
#   "path_tree": "./extdata/Tree/",
#   "tree_names": [
#                 "newseq_LANL_Subtype_01_AE_mafft_aligned_ndrm",
#                 "newseq_LANL_Subtype_A1_mafft_aligned_ndrm"
#                 ],
#   "path_colored_tree":"./extdata/Colored_tree/",
#   "path_clades":"./extdata/Clades/",
#  
#  
#   "local_group": "WA",
#   "focus_group": "KC",
#  
#   "focus_subgroup": null,
#   "index_subgroup": null,
#   "name_subgroup": null,
#  
#   "bs.n": 2,
#   "subtypes": ["01_AE","A1"],
#   "verbose": true
#  }

## ---- echo=TRUE,eval=FALSE,code=Blast_Rcsript-----------------------------------------------------
#  ### Make a blast database ----
#  makeblastdb(LANL_seq,LANL_db,
#              job_index,verbose)
#  # Run the following command for each subtype: (mixing between R and Linux command)
#  #     makeblastdb -in inputs -dbtype nucl -title
#  #       gsub(".fasta","",inputs)
#  #     -out outputs
#  # Where:
#  #     -in: the sequence file to make the database
#  #     -dbtype: the type of the sequences (nucleotide)
#  #     -title : the title of the job
#  #     -out : the name of the output file
#  
#  ### Run blastn ----
#  blastn(Local_seq,LANL_db,max_closest,
#         job_index,verbose)
#  # Run the following command for each subtype: (mixing between R and Linux command)
#  #     blastn -query name_Subtype_new_sequences -db names_db -out
#  #       gsub("_stripped.db","_closest.txt",names_db)
#  #     -max_target_seqs nbre_close -outfmt 6
#  # Where:
#  #     -query: the sequence files ("studied sequences") for which we are looking for the closest in names_db
#  #     -db: the blast database
#  #     -out: the name of the output
#  #     -max_target_seqs: the number of closest sequences to find in the database for each studied sequence

## ---- echo=TRUE,eval=FALSE,code=get_closest_seq_Rcsript-------------------------------------------
#  ### Get closest sequences ----
#  closest_result <- get_closest_seq(LANL_seq,Local_seq,
#                                    path_aligned,seq_names,
#                                    max_closest, # maximal number of closest sequences
#                                    nbre_close,   # number of closest sequences
#                                    min_percent, # minimal percent of closest sequences
#                                    name_HXB2_file,
#                                    job_index,verbose)
#  if(is.na(job_index) || job_index == "NA"){
#   save(closest_result,file=paste(path_aligned,"Closest_result.rda",sep=""))
#  }else{
#   save(closest_result,file=paste(path_aligned,"Closest_result_",job_index,".rda",sep=""))
#  }
#  
#  ### Check that HXB2 is in each file  ----
#  names_file <- paste(path_aligned,seq_names,".fasta",sep="")
#  check_HXB2(names_file,name_HXB2_file,
#             job_index,verbose)

## ---- echo=TRUE,eval=FALSE,code=align_Rcsript-----------------------------------------------------
#  ### Add an outgroup sequences from the closest subtype ----
#  outgroup <- add_outgroup(subtypes,
#                           path_aligned,LANL_seq,seq_names,
#                           outgroup_size,
#                           job_index,verbose)
#  if(is.na(job_index) || job_index == "NA"){
#   save(outgroup,file=paste(path_aligned,"outgroup.rda",sep=""))
#  }else{
#   save(outgroup,file=paste(path_aligned,"outgroup_",job_index,".rda",sep=""))
#  }
#  
#  
#  ### Align ----
#  align_seq(names_file,aligner,thread,
#            job_index,verbose)
#  # Run the following command for each subtype: (mixing between R and Linux command)
#  # mafft --auto --thread 5
#  #      --add tmpdir/names_split[j]
#  #      tmpdir/names_split_aligned[j-1] > tmpdir/names_split_aligned[j]
#  # Where:
#  #     --auto : to get the defaut options in terms of speed and accuracy
#  #     --thread: number of threads (parallele running)
#  #     -add : to merge the output to a given file

## ---- echo=TRUE,eval=FALSE,code=postprocess_align_gaps_Rcsript------------------------------------
#  ### Remove some gaps column ----
#  names_file_aligned <- gsub(x=names_file,".fasta",paste("_",aligner,"_aligned.fasta",sep=""))
#  
#  for(name in names_file_aligned){
#   seq <- read.dna(name,format="fa")
#   seq	<- seq.strip.gap(seq, strip.pc=1-nbre_col_rm_strip/nrow(seq))
#   write.dna(seq,name,format='fasta',colsep='', nbcol=-1)
#  }

## ---- echo=TRUE,eval=FALSE,code=postprocess_align_identical_Rcsript-------------------------------
#  ### Remove the identical sequences ----
#  indexes <- find_identical_sequences(names_file_aligned,path_sequences_meta,
#                                      dna_region,
#                                      job_index,verbose)
#  if(is.na(job_index) || job_index == "NA"){
#   save(indexes,file=paste(path_aligned,"indexes_identical.rda",sep=""))
#  }else{
#   save(indexes,file=paste(path_aligned,"indexes_identical_",job_index,".rda",sep=""))
#  }

## ---- echo=TRUE,eval=FALSE,code=postprocess_align_DRM_Rcsript-------------------------------------
#  ### Remove Drug Resistance Mutations ----
#  rm_DRM(names_file_aligned,name_HXB2,
#         job_index,verbose)

## ---- echo=TRUE,eval=FALSE,code=build_tree_Rcsript------------------------------------------------
#  ### Building trees ----
#  building_tree(names_file_aligned_ndrm,path_tree,
#                tree_builder,
#                bs.n,bs.id,
#                verbose)
#  path_tree <- paste(normalizePath(path_tree),"/",sep="")

## ---- echo=TRUE,eval=FALSE,code=postprocess_tree_Rcsript------------------------------------------
#  ### Postprocess trees (rooting and remove outgroup sequences) ----
#  tree_names <-
#   paste(path_tree,"fasttree_",gsub(x=basename(names_file_aligned_ndrm),".fasta","/"),
#         gsub(x=basename(names_file_aligned_ndrm),".fasta",".newick"),sep='')
#  
#  postprocess_tree(tree_names,path_outgroup,name_HXB2,
#                   bs.n,bs.id,
#                   verbose)

## ---- echo=TRUE,eval=FALSE,code=relabelling_Rcsript-----------------------------------------------
#  ### Relabelling ----
#  tips_relabelling(tree_names,
#                   person,sequences_meta,Country_db,
#                   vars_for_label,
#                   local_group=local_group,
#                   focus_group=focus_group,
#                   focus_subgroup=focus_subgroup,
#                   index_subgroup=index_subgroup,
#                   name_subgroup=name_subgroup,
#                   negative=negative,
#                   local_indicator=local_indicator,
#                   focus_indicator=focus_indicator,
#                   separator=separator,
#                   blacklist_paths=blacklist_paths,
#                   blacklist_names=blacklist_names,
#                   Unknown = Unknown,
#                   local_tip_split=local_tip_split,
#                   external_tip_split=external_tip_split,
#                   local_index_location=local_index_location,
#                   external_index_location=external_index_location,
#                   bs.n=bs.n,bs.id=bs.id,
#                   verbose =verbose)

## ---- echo=TRUE, eval=FALSE-----------------------------------------------------------------------
#   "focus_subgroup": ["HSX_F","HSX_M"], # the categories to focus
#   "index_subgroup": 7, # in which spot is this information in the tip labels
#   "name_subgroup": "HSX" # the name of the subgroup that will be use to relabel the tips

## ---- echo=TRUE,eval=FALSE,code=phyloscanner_Rcsript----------------------------------------------
#  ### Phyloscanner ----
#  wrap_phyloscanner_analyse_tree(renamed_tree_names,path_colored_tree,
#                                 opt_regex=opt_regex,
#                                 bs.n = bs.n,bs.id = bs.id,
#                                 verbose =verbose)

## ---- echo=TRUE,eval=FALSE,code=postprocess_phyloscanner_Rcsript----------------------------------
#  ### Phyloscanner postprocess ----
#  if(!is.null(focus_subgroup)) focus_group <- paste(focus_group,focus_subgroup,sep="")
#  phyloscanner_postprocess(renamed_tree_names,path_colored_tree,
#                           focus_group=focus_group,separator=separator,
#                           bs.n = bs.n,bs.id = bs.id,
#                           verbose =verbose)

## ---- echo=TRUE,eval=FALSE,code=extract_clades_Rcsript--------------------------------------------
#  ### Clades extracting for each bootstrap index ----
#  clades_extracting(renamed_tree_names,
#                    path_colored_tree,path_clades,
#                    focus_group,
#                    bs.n=bs.n,bs.id=bs.id,
#                    verbose=verbose)

