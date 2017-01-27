#' Conversion between gene naming convenctions and organisms
#'
#' @param SourceOrganism Source organism (currently either "human" or "mouse")
#' @param TargetOrganism Target organism (currently either "human" or "mouse")
#' @param Genes Vector of strings containing gene names. NAs can be present and will be removed
#' @param SourceTypes Type of gene names used as input (currently either "Names" or "Ensembl"). Default is "Names".
#' @param TargetTypes Type of gene names to be returned (currently either "Names" or "Ensembl"). Default is "Names". 
#' @param HomologyLevel minimal level of homology (1 by default)
#'
#' @return
#' @export
#'
#' @examples 
#' 
#' GenNames <- c("CCND3", "CD151", "CD2BP2", "CD81", "CDC34", "CDC37", "CDC42BPB", "CDC42EP1","CDC42EP4")
#' 
#' ConvertNames("human", "mouse", SourceTypes = "Names", TargetTypes = "Names", GenNames)
#' 
#' MoN <- ConvertNames("human", "mouse", SourceTypes = "Names", TargetTypes = "Names", GenNames)
#' MoEn <- ConvertNames("human", "mouse", SourceTypes = "Names", TargetTypes = "Ensembl", GenNames)
#' 
#' ConvertNames("mouse", "human", SourceTypes = "Names", TargetTypes = "Names", MoN)
#' ConvertNames("mouse", "human", SourceTypes = "Ensembl", TargetTypes = "Ensembl", MoEn)
#' 
#' 
ConvertNames <- function(SourceOrganism, TargetOrganism, Genes, SourceTypes = "Names", TargetTypes = "Names", HomologyLevel = 1) {

  if(SourceOrganism == "human"){
    
    # Remove NAs
    Genes <- Genes[!is.na(Genes)]
    
    if(TargetOrganism == "mouse"){
      
      if(SourceTypes == "Names"){
        MouseGenes <- biomaRt::getBM(attributes = c("mmusculus_homolog_orthology_confidence",
                                                    "mmusculus_homolog_ensembl_gene",
                                                    "mmusculus_homolog_associated_gene_name"),
                                     filters = "external_gene_name",
                                     values = Genes,
                                     mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
      }
      
      if(SourceTypes == "Ensembl"){
        MouseGenes <- biomaRt::getBM(attributes = c("mmusculus_homolog_orthology_confidence",
                                                    "mmusculus_homolog_ensembl_gene",
                                                    "mmusculus_homolog_associated_gene_name"),
                                     filters = "ensembl_gene_id",
                                     values = Genes,
                                     mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
      }
      
      if(TargetTypes == "Names"){
        return(MouseGenes$mmusculus_homolog_associated_gene_name[MouseGenes$mmusculus_homolog_orthology_confidence >= HomologyLevel &
                                                                   !is.na(MouseGenes$mmusculus_homolog_orthology_confidence)])
      }
      
      if(TargetTypes == "Ensembl"){
        return(MouseGenes$mmusculus_homolog_ensembl_gene[MouseGenes$mmusculus_homolog_orthology_confidence >= HomologyLevel &
                                                                   !is.na(MouseGenes$mmusculus_homolog_orthology_confidence)])
      }
      
    }
    
    if(TargetOrganism == "human"){
      
      if(SourceTypes == "Names"){
        HumanGenes <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                                     filters = "external_gene_name",
                                     values = Genes,
                                     mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
      }
      
      if(SourceTypes == "Ensembl"){
        HumanGenes <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                                     filters = "ensembl_gene_id",
                                     values = Genes,
                                     mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
      }
      
      if(TargetTypes == "Names"){
        return(HumanGenes$external_gene_name[!is.na(HumanGenes$external_gene_name)])
      }
      
      if(TargetTypes == "Ensembl"){
        return(HumanGenes$ensembl_gene_id[!is.na(HumanGenes$ensembl_gene_id)])
      }
      
    }
    
  }
  
  
  if(SourceOrganism == "mouse"){
    
    # Remove NAs
    Genes <- Genes[!is.na(Genes)]
    
    if(TargetOrganism == "human"){
      
      if(SourceTypes == "Names"){
        HumanGenes <- biomaRt::getBM(attributes = c("hsapiens_homolog_orthology_confidence",
                                                    "hsapiens_homolog_ensembl_gene",
                                                    "hsapiens_homolog_associated_gene_name"),
                                     filters = "external_gene_name",
                                     values = Genes,
                                     mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
      }
      
      if(SourceTypes == "Ensembl"){
        HumanGenes <- biomaRt::getBM(attributes = c("hsapiens_homolog_orthology_confidence",
                                                    "hsapiens_homolog_ensembl_gene",
                                                    "hsapiens_homolog_associated_gene_name"),
                                     filters = "ensembl_gene_id",
                                     values = Genes,
                                     mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
      }
      
      if(TargetTypes == "Names"){
        return(HumanGenes$hsapiens_homolog_associated_gene_name[HumanGenes$hsapiens_homolog_orthology_confidence >= HomologyLevel &
                                                                   !is.na(HumanGenes$hsapiens_homolog_orthology_confidence)])
      }
      
      if(TargetTypes == "Ensembl"){
        return(HumanGenes$hsapiens_homolog_ensembl_gene[HumanGenes$hsapiens_homolog_orthology_confidence >= HomologyLevel &
                                                           !is.na(HumanGenes$hsapiens_homolog_orthology_confidence)])
      }
      
    }
    
    
    
    
    if(TargetOrganism == "mouse"){
      
      if(SourceTypes == "Names"){
        MouseGenes <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                                     filters = "external_gene_name",
                                     values = Genes,
                                     mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
      }
      
      if(SourceTypes == "Ensembl"){
        MouseGenes <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                                     filters = "ensembl_gene_id",
                                     values = Genes,
                                     mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
      }
      
      if(TargetTypes == "Names"){
        return(MouseGenes$external_gene_name[!is.na(MouseGenes$external_gene_name)])
      }
      
      if(TargetTypes == "Ensembl"){
        return(MouseGenes$ensembl_gene_id[!is.na(MouseGenes$ensembl_gene_id)])
      }
      
    }
    
  }

}













