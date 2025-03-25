#####################################################################
######################## MGI Validator utils ######################## 
#####################################################################
get_subgraphs <- function(graf, decompose= FALSE){
  # Find nodes with out-degree = 0
  eff_nodes <- V(graf)[degree(graf, mode="out") == 0] %>% names %>%
    gsub(pattern = "_func$", replacement = "", x = .)
  eff_subgraphs <- sapply(eff_nodes, function(eff_node){
    eff_circuit_nodes <- subcomponent(graf, eff_node, mode="in")
    eff_subgraph <- induced_subgraph(graf, eff_circuit_nodes) 
    # if(decompose){
    #   # Find receptor nodes (nodes with in-degree = 0 in the subgraph)
    #   rcp_nodes <- V(graf)[degree(graf, mode="in") == 0]
    #   all_decomposed_subpaths <- lapply(rcp_nodes, function(rcp) {
    #     paths <- all_simple_paths(graf, from = rcp, to = eff_node, mode = "out")
    #     lapply(paths, function(path) induced_subgraph(graf, path)) # Convert paths to subgraphs
    #   })
    #   return(all_decomposed_subpaths)
    # }
    return(eff_subgraph)
  })
  names(eff_subgraphs) <- gsub(pattern = "^N", replacement = "P", x = names(eff_subgraphs)) 
  return(eff_subgraphs)
}

get_not_detectables_from_graph <- function(pathways){
  differences <- lapply(pathways$pathigraphs, function(pathway){
    all_diff <- list()
    mgi_effector.subgraphs <- pathway$effector.subgraphs
    actual_effector.subgraphs <- get_subgraphs(pathway$graph, decompose = F)
    # actual_decomposed.subgraphs <- get_subgraphs(pathway$graph, decompose = T)
    # check if they have same number of subpathways 
    effSub_differences <- union(setdiff(names(actual_effector.subgraphs), names(mgi_effector.subgraphs)),
                                setdiff(names(mgi_effector.subgraphs), names(actual_effector.subgraphs)))
    if(length(effSub_differences)>0){
      warning("differences between actual_effector.subgraphs and mgi_effector subgraphs")
      all_diff$effSub_differences <- effSub_differences
    }
    # Check if the graph has cycles (DAG check)
    # if(!is_dag(pathway$graph))
    #   warning( "Pathway ",pathway$path.id," : ",pathway$path.name,", has cycle(s).")
    subgraphs_vertices <- lapply(actual_effector.subgraphs, get_nodes_from_graph) %>%
      bind_rows %>% distinct 
    graph_vertices <- get_nodes_from_graph(pathway$graph)
    vertices_differences <- setdiff(subgraphs_vertices$name, graph_vertices$name)
    if(length(vertices_differences)>0){
      warning("some vertices will not be detected in subgraph!")
      all_diff$vertices_differences <- vertices_differences
    }
    return(all_diff)
  }) %>% purrr::discard(~ is.null(.) || length(.) == 0)
  return(differences)
}
get_nodes_from_graph <- function(graf){
  return(igraph::as_data_frame(graf, what = "vertice") %>%
           select(c(name, label, genesList)) %>%
           tidyr::unnest(genesList)  %>% 
           filter(genesList!="/" | is.na(genesList)))
}
write_logs <- function(warn_message, log_file = "log.txt", ...){
  write.table(x = paste0(warn_message, "\n##--------------------------------------##"),
              file = log_file, sep = "\t", quote = F, row.names = F, col.names = F, ...)
}

check_integrity_pathway <- function(pathway, log_file = "log.txt", verbose = FALSE) {
  # Get vertices attributes from each graph for a single pathway
  graph_vertices <- get_nodes_from_graph(pathway$graph)
  subgraphs_vertices <- lapply(pathway$subgraphs, get_nodes_from_graph) %>%
    bind_rows %>% distinct
  effector_subgraphs_vertices <- lapply(pathway$effector.subgraphs, get_nodes_from_graph) %>%
    bind_rows %>% distinct
  
  # Initialize list for warnings
  warnings_list <- list()
  
  # Missing nodes in subgraphs/effector subgraphs
  # **Question 1: Are all nodes from the original graph present in the subgraphs?**
  diff_graph_subgraph <- anti_join(graph_vertices, subgraphs_vertices, by = join_by(name, label, genesList))
  if (verbose) cat(sprintf("Checking node integrity in %s: %s...\n", pathway$path.id, pathway$path.name))
  
  if (nrow(diff_graph_subgraph) > 0) {
    warn_message <- sprintf(
      "Pathway %s is missing %d genes and/or %d nodes in subgraphs:\n\t* Genes: %s\n\t* Nodes: %s",
      pathway$path.id, 
      length(unique(diff_graph_subgraph$genesList)), 
      length(unique(diff_graph_subgraph$name)),
      paste(unique(diff_graph_subgraph$genesList), collapse = ", "),
      paste(unique(diff_graph_subgraph$name), collapse = ", ")
    )
    warnings_list <- append(warnings_list, warn_message)
  }
  
  diff_graph_effsubgraph <- anti_join(graph_vertices, effector_subgraphs_vertices, by = join_by(name, label, genesList))
  if (nrow(diff_graph_effsubgraph) > 0) {
    warnings_list <- append(warnings_list, sprintf(
      "Pathway %s is missing %d nodes in effector subgraphs: %s",
      pathway$path.id, nrow(diff_graph_effsubgraph), paste(diff_graph_effsubgraph$name, collapse = ", ")
    ))
  }
  
  # Extra nodes in subgraphs
  # **Question 2: Do subgraphs contain extra nodes not in the original graph?**
  diff_subgraph_graph <- anti_join(subgraphs_vertices, graph_vertices, by = join_by(name, label, genesList))
  if (nrow(diff_subgraph_graph) > 0) {
    warnings_list <- append(warnings_list, sprintf(
      "Subgraphs of pathway %s contain %d extra nodes not in the main graph: %s",
      pathway$path.id, nrow(diff_subgraph_graph), paste(diff_subgraph_graph$name, collapse = ", ")
    ))
  }
  
  diff_effsubgraph_graph <- anti_join(effector_subgraphs_vertices, graph_vertices, by = join_by(name, label, genesList))
  if (nrow(diff_effsubgraph_graph) > 0) {
    warnings_list <- append(warnings_list, sprintf(
      "Effector subgraphs of pathway %s contain %d extra nodes not in the main graph: %s",
      pathway$path.id, nrow(diff_effsubgraph_graph), paste(diff_effsubgraph_graph$name, collapse = ", ")
    ))
  }
  
  #Differences between subgraphs and effector subgraphs
  # **Question 3: Are subgraphs and effector subgraphs identical?**
  diff_subgraph_effsubgraph <- anti_join(subgraphs_vertices, effector_subgraphs_vertices, by = join_by(name, label, genesList))
  if (nrow(diff_subgraph_effsubgraph) > 0) {
    warnings_list <- append(warnings_list, sprintf(
      "Subgraphs of pathway %s contain %d nodes missing in effector subgraphs: %s",
      pathway$path.id, nrow(diff_subgraph_effsubgraph), paste(diff_subgraph_effsubgraph$name, collapse = ", ")
    ))
  }
  
  diff_effsubgraph_subgraph <- anti_join(effector_subgraphs_vertices, subgraphs_vertices, by = join_by(name, label, genesList))
  if (nrow(diff_effsubgraph_subgraph) > 0) {
    warnings_list <- append(warnings_list, sprintf(
      "Effector subgraphs of pathway %s contain %d nodes missing in subgraphs: %s",
      pathway$path.id, nrow(diff_effsubgraph_subgraph), paste(diff_effsubgraph_subgraph$name, collapse = ", ")
    ))
  }
  
  # **Question 4: Are labels in each pathway complete or not?**
  # Checking for discrepancies between `label.id` and the actual labels in the graph
  diff_graph_label.id <- rbind(
    anti_join(graph_vertices[, c("name", "label")], as.data.frame(pathway$label.id), by = join_by(name, label)),
    anti_join(as.data.frame(pathway$label.id), graph_vertices[, c("name", "label")], by = join_by(name, label))
  )
  
  if (nrow(diff_graph_label.id) > 0) {
    warnings_list <- append(warnings_list, sprintf(
      "Labels of pathway %s contain discrepancies between $label.id and the actual labels from the $graph. Affected nodes: %s",
      pathway$path.id, paste(diff_graph_label.id$name, collapse = ", ")
    ))
  }
  
  # Log all warnings if any exist
  if (length(warnings_list) > 0) {
    warning(sprintf("%s : %s has %d warnings!", pathway$path.id, pathway$path.name, length(warnings_list)))
    if(!is.null(log_file )){
      write_logs(warn_message = sprintf("# %s : %s has %d warnings!", pathway$path.id, pathway$path.name, length(warnings_list)),
                 log_file = log_file, append = T)
      write_logs(warn_message = paste("##",warnings_list, collapse = "\n"), 
                 log_file = log_file, append = TRUE)
    }
  }
  ## from MGI iam check if there are some nodes or subpathway that are not detectable easly 
  not_dectable_circuits_and_nodes <- get_not_detectables_from_graph(pathways)
  return(list(
    graph_vertices = graph_vertices,
    subgraphs_vertices = subgraphs_vertices,
    effector_subgraphs_vertices = effector_subgraphs_vertices,
    warnings_list = warnings_list,
    not_dectable_circuits_and_nodes = not_dectable_circuits_and_nodes
  ))
}

check_integrity_MGI <- function(pathways, 
                                log_file = paste0("log_",format(Sys.time(), format= "%Y-%m-%d_%H:%M:%S"),"_.txt"),
                                verbose){
  # log_file <- paste0("log_",format(Sys.time(), format= "%Y-%m-%d_%H:%M:%S"),"_.txt")
  if(!is.null(log_file )) write_logs(warn_message = timestamp(), log_file = log_file, append = F)
  inter_graph_integrity <- lapply(pathways$pathigraphs, check_integrity_pathway,log_file = log_file , verbose= verbose)
  # Here I have to check genes
  all_graph_vertices <- do.call(rbind, lapply(names(inter_graph_integrity), function(pathway_id) {
    inter_graph_integrity[[pathway_id]][["graph_vertices"]]
  })) %>% distinct
  all_subgraphs_vertices <- do.call(rbind, lapply(names(inter_graph_integrity), function(pathway_id) {
    inter_graph_integrity[[pathway_id]][["subgraphs_vertices"]]
  })) %>% distinct
  all_effector_subgraphs_vertices <- do.call(rbind, lapply(names(inter_graph_integrity), function(pathway_id) {
    inter_graph_integrity[[pathway_id]][["effector_subgraphs_vertices"]]
  })) %>% distinct
  all_warnings_list <- lapply(names(inter_graph_integrity), function(pathway_id) {
    inter_graph_integrity[[pathway_id]][["warnings_list"]]
  }) %>% setNames(names(inter_graph_integrity)) %>% purrr::discard(~ is.null(.) || length(.) == 0)
  diff_genes_graph_all.genes <- union(setdiff(pathways$all.genes, all_graph_vertices$genesList), setdiff(all_graph_vertices$genesList, pathways$all.genes)) %>% .[!(. %in% c("NA", NA, "/"))]
  if (length(diff_genes_graph_all.genes) > 0) {
    warn_message <- sprintf(
      "Mismatch detected between pathways$all.genes and the genes present in the graph. \nMissing or extra genes (%d): %s",
      length(diff_genes_graph_all.genes),
      paste(diff_genes_graph_all.genes, collapse = ", ")
    )
    warning(warn_message)
    if(!is.null(log_file )) write_logs(warn_message = warn_message, log_file = log_file, append = T)
    all_warnings_list$diff_genes_graph_all.genes <- warn_message
  }
  diff_genes_subgraph_all.genes <- union(setdiff(pathways$all.genes, all_subgraphs_vertices$genesList), setdiff(all_subgraphs_vertices$genesList, pathways$all.genes)) %>% .[!(. %in% c("NA", NA, "/"))]
  if (length(diff_genes_subgraph_all.genes) > 0) {
    warn_message <- sprintf(
      "Mismatch detected between pathways$all.genes and the genes present in subgraphs. \nMissing or extra genes (%d): %s",
      length(diff_genes_subgraph_all.genes),
      paste(diff_genes_subgraph_all.genes, collapse = ", ")
    )
    warning(warn_message)
    if(!is.null(log_file )) write_logs(warn_message = warn_message, log_file = log_file, append = T)
    all_warnings_list$diff_genes_subgraph_all.genes <- warn_message
  }
  diff_genes_effector_subgraphs_all.genes <- union(setdiff(pathways$all.genes, all_effector_subgraphs_vertices$genesList), setdiff(all_effector_subgraphs_vertices$genesList, pathways$all.genes)) %>% .[!(. %in% c("NA", NA, "/"))]
  if (length(diff_genes_effector_subgraphs_all.genes) > 0) {
    warn_message <- sprintf(
      "Mismatch detected between pathways$all.genes and the genes present in effector subgraphs. \nMissing or extra genes (%d): %s",
      length(diff_genes_effector_subgraphs_all.genes),
      paste(diff_genes_effector_subgraphs_all.genes, collapse = ", ")
    )
    warning(warn_message)
    if(!is.null(log_file )) write_logs(warn_message = warn_message, log_file = log_file, append = T)
    all_warnings_list$diff_genes_effector_subgraphs_all.genes <- warn_message
  }
  return(all_warnings_list)
}
