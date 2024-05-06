library(vidger)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)

# Analysis 1: GO analysis -------------------------------------------------------------
GO_Comparison <- function(Merged_Pvalues, Figure_Path){
  scBSP_Sig_Genes <- Merged_Pvalues$Genes[which(Merged_Pvalues$scBSP<0.05)]
  scBSP_Sig_GO <- enrichGO(scBSP_Sig_Genes,
                           OrgDb = "org.Mm.eg.db", 
                           ont="BP", 
                           keyType="SYMBOL", 
                           pvalueCutoff = 0.0001,
                           readable=TRUE)
  SPARKX_Sig_Genes <- Merged_Pvalues$Genes[which(Merged_Pvalues$SPARKX<0.05)]
  SPARKX_Sig_GO <- enrichGO(SPARKX_Sig_Genes,
                           OrgDb = "org.Mm.eg.db", 
                           ont="BP", 
                           keyType="SYMBOL", 
                           pvalueCutoff = 0.0001,
                           readable=TRUE)
  
  # generate figure
  png(file = Figure_Path,
      width=9, height=12, units = "in", res = 600)
  print(dotplot(scBSP_Sig_GO, showCategory = 20))
  dev.off()
  
  # compare the ranked list of first 100 GO term in scBSP
  
  SPARKX_Rank <- sapply(scBSP_Sig_GO@result$ID[1:100], function(GO_ID){
    if(GO_ID %in% SPARKX_Sig_GO@result$ID){
      return(which(SPARKX_Sig_GO@result$ID == GO_ID))
    }
    else{
      return(NA)
    }
  })
  Output_unit <- c(1 - length(setdiff(scBSP_Sig_GO@result$ID[1:50], SPARKX_Sig_GO@result$ID[1:50])) / 50,
                   1 - length(setdiff(scBSP_Sig_GO@result$ID[1:50], SPARKX_Sig_GO@result$ID[1:100])) / 50,
                   cor.test(1:length(SPARKX_Rank), SPARKX_Rank, method = "kendall")$estimate,
                   cor.test(1:length(SPARKX_Rank), SPARKX_Rank, method = "kendall")$p.value)
  names(Output_unit) <- c("Top50_Diff", "Top50_100_Diff", "Kendall_Tau", "Kendall_Pvalue")
  return(Output_unit)
}



# Analysis 2: Threshold -------------------------------------------------------------

Enriched_Prop <- function(Merged_Pvalues){
  Threshold_Based_Selection <- function(GeneNames, P_Values, Threshold){
    Sig_Genes <- GeneNames[which(P_Values<Threshold)]
    if(length(Sig_Genes) > 0){
      selected_genes_enriched <- enrichGO(Sig_Genes,
                                          OrgDb = "org.Mm.eg.db", 
                                          ont="BP", 
                                          keyType="SYMBOL", 
                                          pvalueCutoff = 0.0001,
                                          readable=TRUE)
      # map enriched GO terms to genes
      if(nrow(as.data.frame(selected_genes_enriched))>0){
        Enriched_Gene_list <- unique(unlist(strsplit(as.data.frame(selected_genes_enriched)$geneID, "/")[[1]]))
        return(c(length(Sig_Genes), length(Enriched_Gene_list)))
      }
      else{
        return(c(length(Sig_Genes), NA))
      }
    }
    else{
      return(c(0, 0))
    }
  }
  
  Rank_Based_Selection <- function(GeneNames, P_Values, Threshold_Rank){
    Sig_Genes <- GeneNames[order(P_Values)[1:Threshold_Rank]]
    selected_genes_enriched <- enrichGO(Sig_Genes,
                                        OrgDb = "org.Mm.eg.db", 
                                        ont="BP", 
                                        keyType="SYMBOL", 
                                        pvalueCutoff = 0.0001,
                                        readable=TRUE)
    # map enriched GO terms to genes
    Enriched_Gene_list <- unique(unlist(strsplit(as.data.frame(selected_genes_enriched)$geneID, "/")[[1]]))
    return(c(length(Sig_Genes), length(Enriched_Gene_list)))
  }
  
  
  
  Results_scBSP_Thres <- sapply(c( 0.005, 0.01, 0.02, 0.03, 0.04, 0.05), function(Threshold_i){
    return(Threshold_Based_Selection(Merged_Pvalues$Genes, Merged_Pvalues$scBSP, Threshold = Threshold_i))
  })
  Results_SPARKX_Thres <- sapply(c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05), function(Threshold_i){
    return(Threshold_Based_Selection(Merged_Pvalues$Genes, Merged_Pvalues$SPARKX, Threshold = Threshold_i))
  })
  
  Results_scBSP_Rank <- sapply(seq(500, 3000, 500), function(Threshold_i){
    return(Rank_Based_Selection(Merged_Pvalues$Genes, Merged_Pvalues$scBSP, Threshold = Threshold_i))
  })
  Results_SPARKX_Rank <- sapply(seq(500, 3000, 500), function(Threshold_i){
    return(Rank_Based_Selection(Merged_Pvalues$Genes, Merged_Pvalues$SPARKX, Threshold = Threshold_i))
  })
  
  
  Results_Thres <- data.frame(Threshold = c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05),
                              scBSP = Results_scBSP_Thres[2,]/Results_scBSP_Thres[1,],
                              SPARKX = Results_SPARKX_Thres[2,]/Results_SPARKX_Thres[1,])
  
  Results_Rank <- data.frame(Threshold = seq(500, 3000, 500),
                             scBSP = Results_scBSP_Rank[2,]/Results_scBSP_Rank[1,],
                             SPARKX = Results_SPARKX_Rank[2,]/Results_SPARKX_Rank[1,])
  
  Output_Unit <- list(Results_Threshold = Results_Thres,
                      Results_Rank = Results_Rank)
  
  return(Output_Unit)
}







# Analysis 3: Threshold baseline -------------------------------------------------------------

Enriched_Prop_Base <- function(Merged_Pvalues){
  Threshold_Based_Selection <- function(GeneNames, P_Values, Threshold){
    Sig_Genes <- GeneNames[which(P_Values<Threshold)]
    if(length(Sig_Genes) > 0){
      Sig_Genes <- sample(GeneNames, length(Sig_Genes))
      selected_genes_enriched <- enrichGO(Sig_Genes,
                                          OrgDb = "org.Mm.eg.db", 
                                          ont="BP", 
                                          keyType="SYMBOL", 
                                          pvalueCutoff = 0.0001,
                                          readable=TRUE)
      # map enriched GO terms to genes
      if(nrow(as.data.frame(selected_genes_enriched))>0){
        Enriched_Gene_list <- unique(unlist(strsplit(as.data.frame(selected_genes_enriched)$geneID, "/")[[1]]))
        return(c(length(Sig_Genes), length(Enriched_Gene_list)))
      }
      else{
        return(c(length(Sig_Genes), NA))
      }
    }
    else{
      return(c(0, 0))
    }
  }
  
  Rank_Based_Selection <- function(GeneNames, P_Values, Threshold_Rank){
    Sig_Genes <- GeneNames[sample(length(GeneNames), Threshold_Rank)]
    selected_genes_enriched <- enrichGO(Sig_Genes,
                                        OrgDb = "org.Mm.eg.db", 
                                        ont="BP", 
                                        keyType="SYMBOL", 
                                        pvalueCutoff = 0.0001,
                                        readable=TRUE)
    # map enriched GO terms to genes
    if(nrow(as.data.frame(selected_genes_enriched))>0){
      Enriched_Gene_list <- unique(unlist(strsplit(as.data.frame(selected_genes_enriched)$geneID, "/")[[1]]))
      return(c(length(Sig_Genes), length(Enriched_Gene_list)))
    }
    else{
      return(c(length(Sig_Genes), NA))
    }
    
    return(c(length(Sig_Genes), length(Enriched_Gene_list)))
  }
  
  
  
  Results_scBSP_Thres <- sapply(c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05), function(Threshold_i){
    return(Threshold_Based_Selection(Merged_Pvalues$Genes, Merged_Pvalues$scBSP, Threshold = Threshold_i))
  })
  Results_SPARKX_Thres <- sapply(c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05), function(Threshold_i){
    return(Threshold_Based_Selection(Merged_Pvalues$Genes, Merged_Pvalues$SPARKX, Threshold = Threshold_i))
  })
  
  Results_scBSP_Rank <- sapply(seq(500, 3000, 500), function(Threshold_i){
    return(Rank_Based_Selection(Merged_Pvalues$Genes, Merged_Pvalues$scBSP, Threshold = Threshold_i))
  })
  Results_SPARKX_Rank <- sapply(seq(500, 3000, 500), function(Threshold_i){
    return(Rank_Based_Selection(Merged_Pvalues$Genes, Merged_Pvalues$SPARKX, Threshold = Threshold_i))
  })
  
  
  Results_Thres <- data.frame(Threshold = c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05),
                              scBSP = Results_scBSP_Thres[2,]/Results_scBSP_Thres[1,],
                              SPARKX = Results_SPARKX_Thres[2,]/Results_SPARKX_Thres[1,])
  
  Results_Rank <- data.frame(Threshold = seq(500, 3000, 500),
                             scBSP = Results_scBSP_Rank[2,]/Results_scBSP_Rank[1,],
                             SPARKX = Results_SPARKX_Rank[2,]/Results_SPARKX_Rank[1,])
  
  Output_Unit <- list(Results_Threshold = Results_Thres,
                      Results_Rank = Results_Rank)
  
  return(Output_Unit)
}

