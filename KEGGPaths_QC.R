#####################################################################
####################### MGI Validator Main Script ###################
#####################################################################
# 
# Description:
# This script validates the integrity of MGI object obtained from KEGG pathways 
# It checks for discrepancies in subgraph, effector subgraph and the main graph
# structures, detects missing or extra nodes, verifies label consistency, 
# and logs potential issues.
#
# Functionality:
# - Extracts effector subgraphs from pathway graphs.
# - Identifies undetectable nodes in subgraphs.
# - Ensures structural consistency between graphs, subgraphs, and effector subgraphs.
# - Logs warnings for missing, extra, or misclassified nodes.
#
# Dependencies:
# - Requires igraph, dplyr, tidyr, and purrr for graph operations.
#
# Output:
# - Warning messages indicating inconsistencies.
# - Log file documenting detected validation issues.
#
# Usage:
# - This script is used within the MGI loaded object from hipathia package 3.0.2 to ensure pathway 
#   models are correctly structured before further analysis.
#
#####################################################################

library("hipathia")
library(igraph)
library(dplyr)
library(purrr)
library(tidyr)
source("src/mgi_validator_functions.R")
pathways <- load_pathways(species = "hsa") 
length(pathways$pathigraphs)
warns <- check_integrity_MGI(pathways, 
                             log_file = paste0("logs/log_",format(Sys.time(), format= "%Y-%m-%d_%H:%M:%S"),"_.txt"),
                             verbose = T)
warns %>% names
# Let's check this example :
##--------------------------------------##
# hsa04020 : Calcium signaling pathway has 2 warnings!
##--------------------------------------##

graph_vertices <- get_nodes_from_graph(pathways$pathigraphs$hsa04020$graph)
subgraphs_vertices <- lapply(pathways$pathigraphs$hsa04020$subgraphs, get_nodes_from_graph) %>%
  bind_rows %>% distinct
missingNodesInSubgraph <- setdiff(graph_vertices$name , subgraphs_vertices) %>% unique
# Check if these node has in/out nodes
missingNodesInSubgraph[degree(pathways$pathigraphs$hsa04020$graph, mode="out")==0] %>% length()
missingNodesInSubgraph[degree(pathways$pathigraphs$hsa04020$graph, mode="in")==0]%>% length()
# check the number of subgraphs 

mgi_effector.subgraphs_hsa04020 <- pathways$pathigraphs$hsa04020$effector.subgraphs
actual_effector.subgraphs_hsa04020 <- get_subgraphs(pathways$pathigraphs$hsa04020$graph, decompose = F)
# names 
mgi_effector.subgraphs_hsa04020 %>% names
actual_effector.subgraphs_hsa04020 %>% names
# nodes
lapply(names(mgi_effector.subgraphs_hsa04020), function(sp){
  all(V(mgi_effector.subgraphs_hsa04020[[sp]])  == V(actual_effector.subgraphs_hsa04020[[sp]]))
  }) %>% unlist() %>% all

# graph 
actual_effector.subgraphs_hsa04020 %>% names %>% length()
"N-hsa04020-6" %in% V(pathways$pathigraphs$hsa04020$graph)$name
lapply(pathways$pathigraphs$hsa04020$effector.subgraphs, function(sp){
  "N-hsa04020-6" %in% V(sp)$name
  }) %>% unlist() %>% any
lapply(actual_effector.subgraphs_hsa04020, function(sp){
  "N-hsa04020-6" %in% V(sp)$name
}) %>% unlist() %>% any
