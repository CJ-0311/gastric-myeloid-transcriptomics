# -*- coding: UTF-8 -*-
# =============================================================================
# PLOS ONE | PONE-D-26-36895 | SINGLE-FILE METHODS CODE
# Integrated single-cell and bulk transcriptomic analysis ...
# Version: 2026-10-25 v2.1 | Windows R | source-time fix
# =============================================================================
# THIS IS A SCIENTIFIC PROVENANCE-AWARE MASTER SCRIPT, NOT A CERTIFICATION THAT
# every newly reconstructed method generated the submitted figures.
#
# SOURCE STATUS:
#   [SAVED]   Reads an existing finalized object/result generated earlier.
#   [REUSED]  Exact pre-existing R scripts embedded as optional legacy functions.
#   [NEW]     Newly written executable implementation of the stated method.
#   [GAP]     Original input/code unavailable; MUST be validated before publication.
#
# HOW TO RUN:
#   source("E:/download/PLOS_ONE_Methods_Fig1_Fig6_AllInOne_v2_1_fixed.R", encoding="UTF-8")
#   methods_status()
#   methods_missing_inputs() # report exact inputs required to finish reproducibility
#   methods_run("M01")   # one module at a time
#   methods_run(c("M01","M02","M04"))
#   # NOTHING is computed automatically on source(); originals are never changed.
#   # After reviewing original legacy code (at end of file), specific optional
#   # functions can be called, but they may contain historical hard-coded paths.
#
# Input conventions for NEW / GAP sections are documented in each function.
# Directory names and cached objects below reflect previous project records.
# =============================================================================

PLOS <- new.env(parent = emptyenv())
PLOS$ROOT <- "E:/download/GSE150290_audit"
PLOS$ATLAS_DIR <- file.path(PLOS$ROOT, "GSE150290_fullcount_integration",
  "Phase5_early_cohort", "Phase7_repaired_and_integrated",
  "Fig1_Original81779_Revised")
PLOS$ATLAS_RDS <- file.path(PLOS$ATLAS_DIR,
  "Original81779_Fig1_Annotated_FINAL.rds")
PLOS$ATLAS_META_CSV <- file.path(PLOS$ATLAS_DIR, "All_Cell_Annotations.csv")
PLOS$MAC_RDS <- file.path(PLOS$ATLAS_DIR, "Fig4_Macrophage_Reanalysis",
  "Phase6_HighResolution", "07_HighRes_Annotated_Candidates.rds")
PLOS$QC_DIR <- file.path(PLOS$ATLAS_DIR, "QC_Methods_Audit_20261021",
  "FourModules_v1")
PLOS$FIG4_DIR <- file.path(PLOS$ROOT, "Fig4_Final", "Fig4_Manuscript_Rebuild_v1")
PLOS$TF_DIR <- file.path(PLOS$FIG4_DIR, "08_TF_SixStates")
PLOS$SIX_STATE_RDS <- file.path(PLOS$TF_DIR,
  "00_TF_Input_1529_DecontX_WithSixStates.rds")
PLOS$REGULON_RDS <- file.path(PLOS$TF_DIR, "01_Reused_SCENIC_Regulons.rds")
PLOS$AUC_RDS <- file.path(PLOS$TF_DIR, "04_SCENIC_Regulon_AUCell_1529.rds")
PLOS$SLINGSHOT_OLD <- file.path(PLOS$ROOT, "Fig4_Final",
  "Fig4C_Trajectory2495_Slingshot.rds")
PLOS$STAGES <- c("NAG", "CAG", "IM", "GC")
PLOS$STAGES_ALL <- c(PLOS$STAGES, "NAT")
PLOS$SEED <- 20261025L
# Output directory is a NEW, never-overwritten folder per run.
PLOS$OUT <- file.path(PLOS$ROOT,
  paste0("E3_Methods_AllInOne_", format(Sys.time(), "%Y%m%d_%H%M%S")))
# For Fig3 differential-expression, this must be a VERIFIED mapping from cells
# or tissue samples to the ORIGINAL 3/3/7/7 analysis units. Do not invent donors.
PLOS$PSEUDOBULK_MAP_CSV <- file.path(PLOS$ROOT, "E3_inputs",
  "Fig3AB_verified_analysis_unit_map.csv")
# For Fig3C (program_name,gene), require the EXACT published gene lists.
PLOS$PROGRAM_GENES_CSV <- file.path(PLOS$ROOT, "E3_inputs",
  "Fig3C_original_program_gene_sets.csv")
# Fig6 input protocol: each cohort needs clinical CSV + sample-by-gene expression CSV.
# Clinical required columns: sample_id, os_months, event (0/1); optional covariates.
# Expression required columns: sample_id + gene symbols; values normalized expression.
PLOS$SURV_DIR <- file.path(PLOS$ROOT, "E3_inputs", "survival")
PLOS$TF_MODULE_CSV <- file.path(PLOS$SURV_DIR, "Fig6_verified_TF_module_targets.csv")
PLOS$COVARIATE_LIST <- list(TCGA=character(), GSE62254=character(),
  GSE84437=character()) # FILL ONLY AFTER confirming original adjusted models.
PLOS$BH_FAMILY_SIZE <- 30L  # Verify against original code; 30 is manuscript claim.

methods_need <- function(pkgs) {
  absent <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(absent)) stop("Install required R packages first: ",
    paste(absent, collapse=", "), call.=FALSE)
}
methods_file <- function(f, meaning="required file") {
  if(!file.exists(f)) stop("Missing ",meaning,": ",f,
    "\nThis is an INPUT/PROVENANCE GAP; no replacement values are invented.",call.=FALSE)
  f
}
methods_dir <- function(tag) {
  p <- file.path(PLOS$OUT,tag)
  dir.create(p,recursive=TRUE,showWarnings=FALSE)
  normalizePath(p,winslash="/",mustWork=TRUE)
}
methods_csv <- function(x,dir,name) {
  stopifnot(is.data.frame(x) || is.matrix(x))
  f <- file.path(dir,name)
  if(file.exists(f)) stop("Refusing to overwrite: ",f)
  utils::write.csv(x,f,row.names=FALSE,na="")
  invisible(f)
}
methods_read <- function(f, required=NULL) {
  z <- utils::read.csv(methods_file(f),check.names=FALSE,
    stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  if(length(required) && !all(required %in% names(z)))
    stop("Missing required columns in ", f, ": ",
      paste(setdiff(required,names(z)),collapse=", "))
  z
}
methods_assay_counts <- function(seu,assay="RNA") {
  methods_need("SeuratObject")
  tryCatch(SeuratObject::GetAssayData(seu,assay=assay,layer="counts"),
    error=function(e) SeuratObject::GetAssayData(seu,assay=assay,slot="counts"))
}
methods_meta <- function(prefer_csv=TRUE) {
  if(prefer_csv && file.exists(PLOS$ATLAS_META_CSV)) {
    d <- methods_read(PLOS$ATLAS_META_CSV)
    if(nrow(d)==81779L && all(c("sample_uid","stage","celltype_fig") %in% names(d))) {
      if(!("cell" %in% names(d))) {
        for(nm in c("cell_barcode","barcode","Cell","X"))
          if(nm %in% names(d)) { d$cell<-as.character(d[[nm]]); break }
      }
      return(d)
    }
  }
  s <- readRDS(methods_file(PLOS$ATLAS_RDS,"final Seurat atlas"))
  stopifnot(ncol(s)==81779L)
  x <- s@meta.data
  x$cell <- rownames(x)
  rm(s);gc(FALSE)
  x
}
methods_load_seu <- function(path=PLOS$ATLAS_RDS) {
  methods_need("SeuratObject")
  s <- readRDS(methods_file(path,"saved Seurat object"))
  if(!inherits(s,"Seurat")) stop("Not a Seurat object: ",path)
  s
}
methods_note <- function(dir,...) {
  cat(paste(...,collapse="\n"),"\n",file=file.path(dir,"PROVENANCE_NOTES.txt"),append=TRUE)
}
methods_capture_env <- function() {
  d<-methods_dir("M00_Environment")
  writeLines(c("E3 methods file:","Result source: previously saved project objects plus new modules",
    "WARNING: New reconstruction is not proof of historical computational execution.",
    capture.output(sessionInfo())),file.path(d,"R_sessionInfo.txt"))
  invisible(d)
}

# ===== M01 | Data sources, QC input and atlas composition [SAVED] =====
methods_M01 <- function() {
  out <- methods_dir("M01_Sources_Preprocessing")
  meta <- methods_meta()
  req <- c("sample_uid","stage","celltype_fig")
  if(!all(req %in% names(meta))) stop("Final atlas metadata is missing: ",
    paste(setdiff(req,names(meta)),collapse=", "))
  if(nrow(meta)!=81779L) stop("Expected 81,779 cells in the final atlas")
  stopifnot(!anyNA(meta[,req]))
  s<-as.data.frame(table(sample_uid=meta$sample_uid,stage=meta$stage),
    stringsAsFactors=FALSE)
  s<-s[s$Freq>0,,drop=FALSE]
  a<-aggregate(Freq~stage,s,sum)
  counts<-aggregate(sample_uid~stage,s,function(x) length(unique(x)))
  names(counts)[2]<-"n_tissue_samples"
  ans<-merge(a,counts,by="stage")
  ord<-match(PLOS$STAGES_ALL,ans$stage)
  ans<-ans[ord[!is.na(ord)],,drop=FALSE]
  methods_csv(s,out,"01_Sample_Cell_Counts.csv")
  methods_csv(ans,out,"02_Stage_Sample_Cell_Counts.csv")
  methods_note(out,
    "Atlas is read from a frozen 81,779-cell Seurat object/metadata CSV.",
    "The original generation of 7 corrected count matrices cannot be reproduced",
    "until their true source/accession/reads provenance and original script are supplied.",
    "GEO anonymous subject metadata != patient identifiers inferred from sample prefixes.")
  print(ans);invisible(ans)
}

# ===== M02 | Seurat/Harmony annotation evidence [SAVED; no fake integration] =====
methods_M02 <- function() {
  out<-methods_dir("M02_Seurat_Harmony_Annotation")
  methods_need("Seurat")
  s<-methods_load_seu()
  stopifnot(ncol(s)==81779L)
  red<-names(s@reductions)
  methods_csv(data.frame(reduction=red),out,"01_Available_Reductions.csv")
  methods_csv(data.frame(metadata_column=names(s@meta.data)),out,
    "02_Metadata_Columns.csv")
  methods_csv(as.data.frame(table(celltype=s$celltype_fig)),out,
    "03_12_Celltype_Counts.csv")
  if("umap" %in% red) {
    methods_need("ggplot2")
    p<-Seurat::DimPlot(s,reduction="umap",group.by="celltype_fig",raster=TRUE) +
      ggplot2::ggtitle("Integrated atlas: 12 major populations")
    ggplot2::ggsave(file.path(out,"Fig1B_Atlas_Annotation_Check.pdf"),p,
      width=10,height=7.5)
  }
  methods_note(out,"Reuses frozen final atlas; DOES NOT pretend to reproduce raw FASTQ",
    "generation, original 2,500-HVG Harmony integration or manual cell-type curation.",
    "Re-run of raw integration requires verified input matrices and batch/sample code.")
  invisible(out)
}

# ===== M03 | Fig1 sample-level 12-cell composition [NEW equivalent summary] ====
methods_M03 <- function() {
  out<-methods_dir("M03_Fig1_Composition")
  meta<-methods_meta()
  req<-c("sample_uid","stage","celltype_fig")
  stopifnot(all(req %in% names(meta)))
  z<-as.data.frame(table(sample_uid=meta$sample_uid,stage=meta$stage,
    celltype=meta$celltype_fig),stringsAsFactors=FALSE)
  # Retain explicit ZERO counts for celltypes absent from an observed sample.
  actual<-unique(meta[,c("sample_uid","stage"),drop=FALSE])
  if(anyDuplicated(actual$sample_uid)) stop("Stage conflicts for tissue sample UID")
  z<-z[paste(z$sample_uid,z$stage) %in%
    paste(actual$sample_uid,actual$stage),,drop=FALSE]
  by_sample<-aggregate(Freq~sample_uid+stage,z,sum)
  z$total<-by_sample$Freq[match(paste(z$sample_uid,z$stage),
    paste(by_sample$sample_uid,by_sample$stage))]
  z$proportion<-z$Freq/z$total
  methods_csv(z,out,"01_Fig1_Sample_Celltype_Proportions.csv")
  methods_need("ggplot2")
  for(type in c("T cells","Macrophages")) {
    w<-z[z$celltype==type,,drop=FALSE]
    if(!nrow(w)) next
    w$stage<-factor(w$stage,levels=PLOS$STAGES_ALL)
    p<-ggplot2::ggplot(w,ggplot2::aes(stage,proportion))+
      ggplot2::geom_boxplot(outlier.shape=NA,width=.5)+
      ggplot2::geom_point(position=ggplot2::position_jitter(width=.07,height=0),size=2)+
      ggplot2::theme_classic(base_size=12)+
      ggplot2::labs(x=NULL,y="Tissue-sample proportion",title=type)
    ggplot2::ggsave(file.path(out,paste0("Fig1_Sample_",gsub(" ","_",type),".pdf")),
      p,width=6.2,height=4.5)
  }
  methods_note(out,"Descriptive tissue-sample proportions, not donor-level inferential tests.",
    "Multiple samples may come from the same individual; no independence is assumed.")
  invisible(z)
}

# ===== M04 | Fig3A/B pseudobulk DE + Hallmark GSEA [NEW; verified map needed] ====
methods_M04 <- function() {
  methods_need(c("SeuratObject","edgeR","limma","scuttle","SummarizedExperiment",
                 "fgsea","msigdbr"))
  out<-methods_dir("M04_Fig3AB_Pseudobulk_GSEA")
  mp<-methods_read(PLOS$PSEUDOBULK_MAP_CSV,c("cell","analysis_unit","stage"))
  stopifnot(!anyDuplicated(mp$cell),!anyNA(mp$analysis_unit),
    all(mp$stage %in% PLOS$STAGES))
  nunit<-unique(mp[,c("analysis_unit","stage")])
  if(anyDuplicated(nunit$analysis_unit)) stop("Pseudobulk unit occurs in multiple stages")
  nums<-table(factor(nunit$stage,levels=PLOS$STAGES))
  if(!identical(as.integer(nums),c(3L,3L,7L,7L)))
    stop("Fig3AB must use original verified 3/3/7/7 analysis units. Found: ",
      paste(nums,collapse="/"))
  s<-methods_load_seu()
  ct<-methods_assay_counts(s)
  if(!all(mp$cell %in% colnames(ct)))stop("Analysis map contains unknown cell barcodes")
  cts<-ct[,mp$cell,drop=FALSE]
  groups<-factor(mp$analysis_unit,levels=unique(nunit$analysis_unit))
  ag<-scuttle::aggregateAcrossCells(cts,ids=groups,statistics="sum")
  mat<-if(inherits(ag,"SummarizedExperiment")) {
    SummarizedExperiment::assay(ag,"sum")
  } else as.matrix(ag)
  if(is.null(colnames(mat)) || !setequal(colnames(mat),levels(groups)))
    colnames(mat)<-levels(groups)
  unit_st<-nunit$stage[match(colnames(mat),nunit$analysis_unit)]
  if(anyNA(unit_st))stop("Pseudobulk columns failed to map to stage")
  y<-edgeR::DGEList(counts=mat,group=unit_st)
  design<-model.matrix(~0+factor(unit_st,levels=PLOS$STAGES))
  colnames(design)<-PLOS$STAGES
  keep<-edgeR::filterByExpr(y,design=design)
  y<-y[keep,,keep.lib.sizes=FALSE]
  y<-edgeR::calcNormFactors(y,method="TMM")
  v<-limma::voom(y,design,plot=FALSE)
  fit<-limma::lmFit(v,design)
  contr<-limma::makeContrasts(
    NAG=NAG-(CAG+IM+GC)/3,CAG=CAG-(NAG+IM+GC)/3,
    IM=IM-(NAG+CAG+GC)/3,GC=GC-(NAG+CAG+IM)/3,levels=design)
  fit<-limma::eBayes(limma::contrasts.fit(fit,contr))
  methods_csv(nunit,out,"01_Analysis_Units.csv")
  for(st in PLOS$STAGES) {
    tt<-limma::topTable(fit,coef=st,number=Inf,sort.by="none")
    tt$gene<-rownames(tt)
    methods_csv(tt,out,paste0("Fig3A_DE_",st,".csv"))
    ranking<-tt$t;names(ranking)<-tt$gene
    ranking<-sort(ranking[is.finite(ranking) & !duplicated(names(ranking))],decreasing=TRUE)
    h<-tryCatch(msigdbr::msigdbr(species="Homo sapiens",collection="H"),
      error=function(e)msigdbr::msigdbr(species="Homo sapiens",category="H"))
    terms<-split(h$gene_symbol,h$gs_name)
    fg<-fgsea::fgseaMultilevel(pathways=terms,stats=ranking,minSize=10,maxSize=500)
    df<-as.data.frame(fg)
    if("leadingEdge" %in% names(df))
      df$leadingEdge<-vapply(df$leadingEdge,paste,character(1),collapse=";")
    methods_csv(df,out,paste0("Fig3B_Hallmark_GSEA_",st,".csv"))
  }
  methods_note(out,"NEW reconstructed analysis, run ONLY on verified original unit mapping.",
    "Non-independent pseudobulk units render DE P/FDR exploratory; no donor-level inference.",
    "Requires original feature-level filters, contrasts and gene-set version confirmation.")
  invisible(out)
}

# ===== M04C | Fig3C program scores at 3/3/7/9 tissue samples [NEW] =====
methods_M04C <- function() {
  methods_need(c("SeuratObject","scuttle","SummarizedExperiment","edgeR","ggplot2"))
  out<-methods_dir("M04C_Fig3C_Sample_Programs")
  gs<-methods_read(PLOS$PROGRAM_GENES_CSV,c("program","gene"))
  s<-methods_load_seu()
  m<-s@meta.data
  stopifnot(all(c("sample_uid","stage") %in% names(m)))
  use<-which(as.character(m$stage) %in% PLOS$STAGES)
  sample_tag<-as.character(m$sample_uid[use]);smeta<-unique(data.frame(
    sample_uid=sample_tag,stage=as.character(m$stage[use])))
  if(anyDuplicated(smeta$sample_uid))stop("Sample UID has multiple disease stages")
  cnt<-table(factor(smeta$stage,levels=PLOS$STAGES))
  if(!identical(as.integer(cnt),c(3L,3L,7L,9L)))
    stop("Fig3C requires 3/3/7/9 tissue samples; found ",paste(cnt,collapse="/"))
  ct<-methods_assay_counts(s)[,use,drop=FALSE]
  u<-factor(sample_tag,levels=smeta$sample_uid)
  agg<-scuttle::aggregateAcrossCells(ct,ids=u,statistics="sum")
  mat<-if(inherits(agg,"SummarizedExperiment")) {
    SummarizedExperiment::assay(agg,"sum")
  } else as.matrix(agg)
  colnames(mat)<-smeta$sample_uid
  d<-edgeR::DGEList(counts=mat)
  d<-edgeR::calcNormFactors(d,method="TMM")
  logcpm<-edgeR::cpm(d,log=TRUE,prior.count=1)
  common<-unique(gs$gene[gs$gene %in% rownames(logcpm)])
  if(!length(common))stop("No program genes matched expression data")
  z<-t(scale(t(logcpm[common,,drop=FALSE])))
  # A gene with zero across-sample variance cannot have a finite Z score.
  dropped<-rownames(z)[rowSums(!is.finite(z))>0]
  methods_csv(data.frame(gene=dropped),out,"00_Program_ZeroVariance_Genes_Excluded.csv")
  z<-z[rowSums(!is.finite(z))==0,,drop=FALSE]
  if(!nrow(z))stop("No program genes with finite across-sample Z-scores")
  # Explicit names are derived from the gene-set keys (not from a guessed title).
  programs<-split(gs$gene,gs$program)
  vals<-lapply(names(programs),function(k){
    genes<-intersect(unique(programs[[k]]),rownames(z))
    if(!length(genes))return(NULL)
    data.frame(sample_uid=colnames(z),program=k,
      score=colMeans(z[genes,,drop=FALSE]),n_genes=length(genes))
  })
  res<-do.call(rbind,vals)
  if(is.null(res))stop("None of the verified program gene sets matched")
  res$stage<-smeta$stage[match(res$sample_uid,smeta$sample_uid)]
  methods_csv(res,out,"Fig3C_TissueSample_ProgramScores.csv")
  p<-ggplot2::ggplot(res,ggplot2::aes(stage,score))+
    ggplot2::geom_boxplot(outlier.shape=NA,width=.5)+
    ggplot2::geom_point(position=ggplot2::position_jitter(width=.08,height=0),size=1.9)+
    ggplot2::facet_wrap(~program,scales="free_y")+
    ggplot2::theme_classic(base_size=12)+
    ggplot2::labs(y="Mean gene-wise Z-score",x=NULL)
  ggplot2::ggsave(file.path(out,"Fig3C_Program_Scores_Check.pdf"),p,width=11,height=4.5)
  methods_note(out,"NEW reconstruction using VERIFIED program-gene list and tissue-sample counts.",
    "Fig3C has 3/3/7/9 tissue sample points; distinct from Fig3AB 3/3/7/7 units.")
  invisible(res)
}

# ===== M05 | CellChat saved-stage networks [REUSED computational definition] ==
methods_M05 <- function() {
  methods_need("CellChat")
  out<-methods_dir("M05_CellChat_Saved_Audit")
  ccdir<-file.path(PLOS$ATLAS_DIR,"Fig2_CellChat_Revised")
  paths<-file.path(ccdir,paste0("CellChat_",PLOS$STAGES,"_primary.rds"))
  if(!all(file.exists(paths)))stop("Missing CellChat saved objects: ",
    paste(paths[!file.exists(paths)],collapse=";"))
  rows<-list(); params<-list()
  for(i in seq_along(paths)){
    st<-PLOS$STAGES[i];o<-readRDS(paths[i]);w<-as.matrix(o@net$weight)
    if(!identical(nrow(w),ncol(w)))stop("Non-square CellChat weight matrix")
    if(!identical(rownames(w),colnames(w)))
      w<-w[,rownames(w),drop=FALSE]
    vout<-rowSums(w);vin<-colSums(w)
    rows[[st]]<-data.frame(stage=st,celltype=rownames(w),
      incoming_weight=unname(vin),outgoing_weight=unname(vout),
      incoming_rank=rank(-vin,ties.method="min"),
      outgoing_rank=rank(-vout,ties.method="min"),
      n_celltypes=nrow(w))
    pars<-o@options$parameter
    getpar<-function(nm)if(is.null(pars[[nm]]))NA_character_ else
      paste(as.character(pars[[nm]]),collapse=";")
    params[[st]]<-data.frame(stage=st,n_cells=length(o@idents),
      type.mean=getpar("type.mean"),raw.use=getpar("raw.use"),
      population.size=getpar("population.size"),nboot=getpar("nboot"),
      Kh=getpar("Kh"),min.cells="UNVERIFIED",DB_version="UNVERIFIED")
  }
  methods_csv(do.call(rbind,rows),out,"Fig2_WithinStage_CellChat_Centrality.csv")
  methods_csv(do.call(rbind,params),out,"CellChat_Saved_Settings.csv")
  methods_note(out,"Weighted network strengths are within-stage descriptive quantities.",
    "Stage networks have different celltype counts/cell numbers; no inferential stage test.",
    "Original min.cells and precise CellChatDB subset still need upstream scripts.")
  invisible(out)
}

# ===== M06 | Macrophage reclustering/marker inventory [NEW re-run] ============
methods_M06 <- function() {
  methods_need(c("Seurat","ggplot2"))
  out<-methods_dir("M06_Macrophage_Seurat_Rerun")
  x<-methods_load_seu(PLOS$MAC_RDS)
  if(ncol(x)!=3101L)stop("Unexpected number of macrophage candidate cells")
  Seurat::DefaultAssay(x)<-"RNA"
  # This is an independently regenerated 17-cluster verification, NOT silently
  # the provenance of the finalized six groups or Fig4A-E.
  x<-Seurat::NormalizeData(x,verbose=FALSE)
  x<-Seurat::FindVariableFeatures(x,nfeatures=2500,verbose=FALSE)
  v<-Seurat::VariableFeatures(x)
  x<-SeuratObject::`VariableFeatures<-`(x,value=v[!grepl("^(MT-|RPL|RPS)",v)])
  x<-Seurat::ScaleData(x,verbose=FALSE)
  x<-Seurat::RunPCA(x,verbose=FALSE)
  x<-Seurat::FindNeighbors(x,dims=1:20,verbose=FALSE)
  x<-Seurat::FindClusters(x,resolution=1.0,verbose=FALSE)
  x<-Seurat::RunUMAP(x,dims=1:20,seed.use=PLOS$SEED,verbose=FALSE)
  mark<-Seurat::FindAllMarkers(x,only.pos=FALSE,min.pct=.20,logfc.threshold=.25)
  methods_csv(mark,out,"Independent_3101_HighRes_Markers.csv")
  saveRDS(x,file.path(out,"Independent_3101_HighRes_Rerun.rds"))
  methods_note(out,"NEW rerun, not a claim that cell labels exactly match the original 17 clusters.",
    "Final 1529-cell six-state annotations must be sourced from the saved audited input.")
  invisible(out)
}

# ===== M07 | Doublets, ambient RNA, robustness [SAVED audit] ==================
methods_M07 <- function() {
  out<-methods_dir("M07_scDblFinder_DecontX_QC")
  files<-c(
    mapping=file.path(PLOS$QC_DIR,"03_Stability","01_All3101_Audit_Mapping.csv"),
    retention=file.path(PLOS$QC_DIR,"03_Stability","02_State_Retention_Sensitivity.csv"),
    marker=file.path(PLOS$QC_DIR,"03_Stability","04_State_Marker_Retention.csv"))
  for(f in files)methods_file(f,"QC saved audit")
  m<-methods_read(files[["mapping"]],c("cell","cluster","doublet_class","contamination"))
  stopifnot(nrow(m)==3101L,!anyDuplicated(m$cell))
  s<-data.frame(scenario=c("Original3101","Singlets","Singlets_and_DecontX_le_0.2"),
    n_cells=c(nrow(m),sum(m$doublet_class!="doublet"),
      sum(m$doublet_class!="doublet" & m$contamination<=.2)))
  methods_csv(s,out,"01_3101_QC_Retention.csv")
  a<-methods_read(files[["retention"]]);methods_csv(a,out,"02_FiveState_Sensitivity.csv")
  methods_csv(methods_read(files[["marker"]]),out,"03_Original_Marker_Sensitivity.csv")
  if(file.exists(file.path(PLOS$QC_DIR,"03_Stability","07_AdjustedRandIndex_ClusterStability.csv")))
    methods_csv(methods_read(file.path(PLOS$QC_DIR,"03_Stability",
      "07_AdjustedRandIndex_ClusterStability.csv")),out,"04_ARI.csv")
  methods_note(out,"READS precomputed scDblFinder/DecontX/QC summaries; no re-fitting.",
    "DecontX marker detection is assessed using >=0.5 corrected UMI threshold.",
    "The originally verified 5-state QC GO terms are not proof for the mixed-marker 6th state.")
  invisible(out)
}

# ===== M08 | CIBERSORT saved bulk/scRNA deconvolution audit [SAVED] ============
methods_M08 <- function() {
  out<-methods_dir("M08_CIBERSORT_Saved_Results")
  base<-file.path(PLOS$ATLAS_DIR,"Fig3_ABC_Revised")
  paths<-c(
    bulk=file.path(base,"Fig3D_Bulk_Normalized10Celltypes.csv"),
    trends=file.path(base,"Fig3E_CrossPlatform_Trends.csv"),
    samples=file.path(base,"Fig3E_CrossPlatform_Sample_Data.csv"))
  found<-paths[file.exists(paths)]
  if(!length(found))stop("No saved CIBERSORT results found: ",base)
  for(nm in names(found))methods_csv(methods_read(found[[nm]]),out,
    paste0("SAVED_",nm,".csv"))
  methods_note(out,"This module audits pre-existing CIBERSORT outputs; does NOT recompute",
    "deconvolution without the original reference signature and CIBERSORT code.",
    "Four immune lineages are re-normalized within their selected denominator;",
    "this is different from a proportion of all single cells.")
  invisible(out)
}

# ===== M09 | FCN1 root plus alternative Slingshot fits [NEW sensitivity] =======
methods_M09 <- function() {
  methods_need(c("slingshot","SeuratObject","ggplot2"))
  out<-methods_dir("M09_Slingshot_Root_Sensitivity")
  x<-methods_load_seu(PLOS$SIX_STATE_RDS)
  stopifnot(ncol(x)==1529L)
  meta<-x@meta.data
  col_state<-intersect(c("fig4_state_display","Fig4_six_state"),names(meta))[1]
  if(is.na(col_state))stop("Saved 1529-cell object has no six-state annotation")
  states<-as.character(meta[[col_state]])
  if(!all(!is.na(states)))stop("Missing six-state labels")
  if(!("umap" %in% names(x@reductions)))
    stop("No saved UMAP in 1529-cell input. Do not substitute a different embedding")
  xy<-SeuratObject::Embeddings(x[["umap"]])[,1:2,drop=FALSE]
  mixed<-unique(states[grepl("HLA-II|GAST",states)])
  if(length(mixed)!=1L)stop("Ambiguous mixed-marker state identifier")
  fc<-unique(states[grepl("FCN1",states)])
  if(length(fc)!=1L)stop("Cannot resolve FCN1 starting state")
  configs<-list(FCN1_six=list(cells=rownames(xy),root=fc),
                HLA_GAST_six=list(cells=rownames(xy),root=mixed),
                FCN1_five=list(cells=rownames(xy)[states!=mixed],root=fc))
  z<-list()
  for(nm in names(configs)) {
    spec<-configs[[nm]];keep<-spec$cells
    lab<-states[match(keep,rownames(xy))]
    fit<-slingshot::slingshot(xy[keep,,drop=FALSE],clusterLabels=lab,
      start.clus=spec$root)
    saveRDS(fit,file.path(out,paste0("Slingshot_",nm,".rds")))
    ln<-slingshot::slingLineages(fit)
    z[[nm]]<-data.frame(scenario=nm,root=spec$root,
      lineage=seq_along(ln),path=vapply(ln,paste,character(1),collapse=" -> "))
  }
  methods_csv(do.call(rbind,z),out,"Root_Sensitivity_Lineage_Paths.csv")
  methods_note(out,"These are NEW starting-point sensitivity fits, not the historical Fig4D object.",
    "Original 3101-cell vs 1529-cell UMAP coordinates are NOT interchangeable.",
    "No biological developmental direction is inferred from stage enrichment alone.")
  invisible(out)
}

# ===== M10 | SCENIC 326 regulons AUCell and final six-state summaries [SAVED] ==
methods_M10 <- function() {
  methods_need(c("SeuratObject","ggplot2"))
  out<-methods_dir("M10_SCENIC_AUCell_Regulons")
  x<-methods_load_seu(PLOS$SIX_STATE_RDS)
  auc<-readRDS(methods_file(PLOS$AUC_RDS,"326-regulon AUCell cache"))
  if(inherits(auc,"aucellResults")) {
    methods_need("AUCell");auc<-AUCell::getAUC(auc)
  }
  auc<-as.matrix(auc)
  if(!setequal(colnames(auc),colnames(x)))stop("AUCell cells != six-state cells")
  auc<-auc[,colnames(x),drop=FALSE]
  if(ncol(auc)!=1529L || nrow(auc)!=326L)
    stop("Unexpected regulon/cell dimensions: ",nrow(auc),"/",ncol(auc))
  m<-x@meta.data
  id<-intersect(c("fig4_state_display","Fig4_six_state"),names(m))[1]
  if(is.na(id))stop("No six-state labels in input Seurat object")
  state<-as.character(m[[id]])
  tf_names<-c("REL","MAF","STAT1","CEBPB","FOSL2","IRF1","NFE2L2",
              "JUNB","SPI1","PRDM1","BHLHE40")
  info<-data.frame(regulon=rownames(auc),
    TF=sub("(_extended)? \\([0-9]+g\\)$","",rownames(auc)),
    extended=grepl("_extended \\(",rownames(auc)))
  info$n_targets<-as.integer(sub("^.*\\(([0-9]+)g\\)$","\\1",info$regulon))
  choose<-lapply(tf_names,function(tf){
    d<-info[info$TF==tf & is.finite(info$n_targets) &
      info$n_targets>=20 & info$n_targets<=1000,,drop=FALSE]
    if(!nrow(d))return(NULL)
    d[order(d$extended,d$n_targets),,drop=FALSE][1,,drop=FALSE]
  })
  choose<-do.call(rbind,choose)
  if(!nrow(choose))stop("None of the prespecified regulons were found")
  methods_csv(choose,out,"Regulon_Selection_Audit.csv")
  X<-t(auc[choose$regulon,,drop=FALSE]);colnames(X)<-choose$TF
  state_means<-rowsum(X,state)/as.numeric(table(state)[rownames(rowsum(X,state))])
  long<-data.frame(expand.grid(state=rownames(state_means),TF=colnames(state_means)))
  long$mean_auc<-as.vector(state_means)
  long$within_regulon_z<-ave(long$mean_auc,long$TF,FUN=function(v){
    s<-sd(v);if(!is.finite(s)||s<1e-10)return(rep(0,length(v)))
    pmax(-2,pmin(2,(v-mean(v))/s))})
  methods_csv(long,out,"Fig5A_SixState_AUCell_Activity.csv")
  if(!all(c("sample_uid","stage") %in% names(m)))
    stop("Saved cell metadata lacks sample_uid or stage")
  # Fig5C state-centering: Z-score each regulon across cells, subtract
  # six-state mean and average by existing SAMPLE UNIT. This is NOT an
  # independent-patient test. Hist GC collapse requires verified mapping.
  TF<-intersect(c("REL","MAF","STAT1","CEBPB"),colnames(X))
  Xz<-apply(X[,TF,drop=FALSE],2L,function(v){
    s<-sd(v);if(!is.finite(s)||s<1e-10)return(rep(0,length(v)))
    (v-mean(v))/s})
  Xz<-as.matrix(Xz);rownames(Xz)<-rownames(X)
  zmeans<-rowsum(Xz,state)/as.numeric(table(state)[rownames(rowsum(Xz,state))])
  Xr<-Xz-zmeans[state,,drop=FALSE]
  smp<-as.character(m$sample_uid)
  by_unit<-rowsum(Xr,smp)/as.numeric(table(smp)[rownames(rowsum(Xr,smp))])
  u<-data.frame(sample_uid=rownames(by_unit),n_cells=as.numeric(table(smp)[rownames(by_unit)]),
    stage=as.character(m$stage[match(rownames(by_unit),smp)]))
  res<-do.call(rbind,lapply(seq_len(ncol(by_unit)),function(j){
    data.frame(u,TF=colnames(by_unit)[j],state_adjusted_Z=by_unit[,j])}))
  res<-res[res$n_cells>=5L,,drop=FALSE]
  methods_csv(res,out,"Fig5C_SampleUnit_StateAdjusted_AUCell.csv")
  methods_note(out,"SCENIC regulons inferred previously with hg19 motif database; AUCell reused.",
    "Fig5A summary and Fig5C SAMPLE-unit audit are reconstructed from cached 326x1529 data.",
    "Original GC patient-level grouping and exact final Fig5D/E community rules require confirmation.",
    "Target number and AUCell score are distinct variables.")
  invisible(out)
}

# ===== M11 | Fig5D/E regulator modules + GO outputs [GAP-aware] ===============
methods_M11 <- function() {
  out<-methods_dir("M11_TF_Modules_GO")
  rds<-PLOS$REGULON_RDS
  reg<-readRDS(methods_file(rds,"saved original SCENIC regulon list"))
  if(!is.list(reg))stop("SCENIC regulon input must be a named list")
  methods_csv(data.frame(regulon=names(reg),target_count=lengths(reg)),
    out,"Regulon_Target_Counts.csv")
  # Published selected-module inputs are required. Do NOT choose a community
  # algorithm or 3 modules post hoc and call it the original implementation.
  f<-file.path(PLOS$TF_DIR,"Fig5_Final_Modules_Verified.csv")
  if(!file.exists(f)){
    methods_note(out,"MISSING verified final Fig5D/E module membership file.",
      "Required columns: TF, module; also original network edges/thresholds/GO inputs.",
      "Historical full and top-100 graphs had differing community assignments (12/13).")
    message("Fig5D/E provenance GAP recorded. No fictitious community algorithm used.")
    return(invisible(out))
  }
  z<-methods_read(f,c("TF","module"))
  if(anyDuplicated(z$TF))stop("Ambiguous TF-to-module mapping")
  methods_csv(z,out,"Verified_TF_Module_Membership.csv")
  methods_note(out,"Verified module membership imported. GO code requires exact target-union background",
    "and module-specific gene universe from the original enrichment script.")
  invisible(out)
}

# ===== M12 | Fig6 discovery Cox + 2 external cohorts + KM [NEW; raw input gap] =
methods_M12 <- function() {
  methods_need("survival")
  out<-methods_dir("M12_TCGA_External_Cox")
  tf<-methods_read(PLOS$TF_MODULE_CSV,c("module","gene"))
  # This input is frozen, verified MODULE target-gene table (26 TF / 4387 targets),
  # NOT constructed by guessing from incomplete 20/26 TF AUCell coverage.
  cohorts<-c("TCGA","GSE62254","GSE84437")
  pair<-function(st){
    cl<-methods_read(file.path(PLOS$SURV_DIR,paste0(st,"_clinical.csv")),
      c("sample_id","os_months","event"))
    ex<-methods_read(file.path(PLOS$SURV_DIR,paste0(st,"_expression.csv")),"sample_id")
    if(anyDuplicated(cl$sample_id)||anyDuplicated(ex$sample_id))
      stop("Duplicate sample_id; original de-duplication procedure must be verified for ",st)
    z<-merge(cl,ex,by="sample_id")
    z<-z[is.finite(z$os_months) & z$os_months>0 & z$event %in% c(0,1),,drop=FALSE]
    if(!nrow(z))stop("No matching eligible survival samples in ",st)
    z
  }
  dat<-setNames(lapply(cohorts,pair),cohorts)
  counts<-data.frame(cohort=cohorts,n=sapply(dat,nrow),
    events=sapply(dat,function(z)sum(z$event)))
  methods_csv(counts,out,"01_Eligible_Cohort_N_and_Events.csv")
  # Discovery gene-level Cox; keep tests distinct from clinical multivariable.
  cox1<-function(z,gene,vars=character()){
    if(!(gene %in% names(z)))return(NULL)
    v<-suppressWarnings(as.numeric(z[[gene]]))
    if(sum(is.finite(v))<10||length(unique(v[is.finite(v)]))<2)return(NULL)
    d<-data.frame(os_months=z$os_months,event=z$event,x=v)
    covs<-intersect(vars,names(z))
    for(cc in covs)d[[cc]]<-z[[cc]]
    d<-d[stats::complete.cases(d),,drop=FALSE]
    if(nrow(d)<10 || sum(d$event)<4)return(NULL)
    formula<-stats::as.formula(paste0("survival::Surv(os_months,event) ~ x",
      if(length(covs))paste0(" + ",paste(covs,collapse=" + ")) else ""))
    fit<-tryCatch(survival::coxph(formula,data=d),error=function(e)NULL)
    if(is.null(fit))return(NULL)
    a<-summary(fit)$coefficients
    ci<-summary(fit)$conf.int
    if(!("x" %in% rownames(a)))return(NULL)
    data.frame(gene=gene,n=nrow(d),events=sum(d$event),
      HR=unname(ci["x","exp(coef)"]),
      HR_lo=unname(ci["x","lower .95"]),HR_hi=unname(ci["x","upper .95"]),
      P=unname(a["x","Pr(>|z|)"]))
  }
  all_genes<-sort(unique(tf$gene))
  discovery<-do.call(rbind,lapply(all_genes,function(g)cox1(dat$TCGA,g)))
  if(is.null(discovery))stop("No testable module targets in discovery expression data")
  # This NEW prioritization uses Cox P rank; compare with ORIGINAL 30-gene list
  # before claiming to reproduce Fig6A. Freeze the 30 candidates independently.
  tab<-merge(unique(tf[,c("module","gene")]),discovery,by="gene")
  fixed<-do.call(rbind,lapply(split(tab,tab$module),function(z){
    z<-z[order(z$P,z$gene),,drop=FALSE]
    head(z,10)}))
  if(nrow(fixed)!=30L)stop("Expected 10 candidates from each of 3 modules; found ",nrow(fixed))
  if(anyDuplicated(fixed$gene)) stop("Repeated genes across modules: validate original deduplication policy")
  methods_csv(discovery,out,"02_Discovery_All_Testable_Cox.csv")
  methods_csv(fixed,out,"03_New_30_Candidate_PValue_Priority.csv")
  results<-list()
  for(st in cohorts){
    d<-dat[[st]]
    tab<-do.call(rbind,lapply(fixed$gene,function(g)cox1(d,g)))
    if(is.null(tab))next
    tab$cohort<-st
    # Fixed n=30 controls the P value family ONLY if original used this.
    # R's p.adjust supports n>length(p) and treats missing tests as p=1.
    tab$BH_fixed30<-stats::p.adjust(tab$P,method="BH",n=PLOS$BH_FAMILY_SIZE)
    tab$BH_evaluable<-stats::p.adjust(tab$P,method="BH")
    results[[st]]<-tab
  }
  r<-do.call(rbind,results)
  methods_csv(r,out,"04_External_And_TCGA_Cox.csv")
  # Adjusted analysis is executed ONLY when verified clinical variable list given.
  adjusted<-list()
  for(st in cohorts){
    cv<-PLOS$COVARIATE_LIST[[st]]
    if(!length(cv))next
    if(!all(cv %in% names(dat[[st]])))stop("Missing adjustment variables in ",st)
    a<-do.call(rbind,lapply(fixed$gene,function(g)cox1(dat[[st]],g,cv)))
    if(!is.null(a)){a$cohort<-st;adjusted[[st]]<-a}
  }
  if(length(adjusted))methods_csv(do.call(rbind,adjusted),out,
    "05_Adjusted_Cox_Verified_Covariates.csv")
  # Descriptive independent consistency check, not an externalized risk model.
  ok<-list()
  for(g in unique(fixed$gene)) {
    x<-r[r$gene==g,,drop=FALSE]
    if(nrow(x)!=3 || !setequal(x$cohort,cohorts))next
    if(all(x$P<0.05) && length(unique(sign(log(x$HR))))==1L)
      ok[[g]]<-data.frame(gene=g,direction=ifelse(all(x$HR>1),"adverse","protective"),
        external_bh_fixed30=all(x$BH_fixed30[x$cohort!="TCGA"]<.05))
  }
  if(length(ok)) methods_csv(do.call(rbind,ok),out,
    "06_Exploratory_CrossCohort_Gene_Consistency.csv")
  # Saved per-cohort KM curves for the SAME three representative genes.
  # New calculation for auditing only; NOT a claim of identical original KM.
  reps<-intersect(c("IGFBP7","SPARC","IKBKE"),fixed$gene)
  if(length(reps)){
    pdf(file.path(out,"07_Three_Representative_Genes_KM_Audit.pdf"),
      width=11,height=11,onefile=TRUE)
    oldpar<-par(mfrow=c(3,3),mar=c(4,4,3,1))
    kmout<-list()
    for(g in c("IGFBP7","SPARC","IKBKE"))for(st in cohorts){
      z<-dat[[st]]
      if(!(g %in% names(z))){plot.new();title(main=paste(st,g,"not measured"));next}
      val<-suppressWarnings(as.numeric(z[[g]]))
      v<-is.finite(val)&is.finite(z$os_months)&z$os_months>0&z$event %in% c(0,1)
      if(sum(v)<10){plot.new();title(main=paste(st,g,"insufficient"));next}
      zz<-data.frame(time=z$os_months[v],event=z$event[v],expression=val[v])
      med<-median(zz$expression)
      zz$group<-factor(ifelse(zz$expression>med,"High","Low"),levels=c("Low","High"))
      if(length(unique(zz$group))!=2){plot.new();title(main=paste(st,g,"no split"));next}
      fit<-survival::survfit(survival::Surv(time,event)~group,data=zz)
      lr<-survival::survdiff(survival::Surv(time,event)~group,data=zz)
      pp<-stats::pchisq(lr$chisq,df=1,lower.tail=FALSE)
      plot(fit,col=c("#4EAB8D","#C55286"),lwd=2,mark.time=FALSE,
        xlab="Overall survival (months)",ylab="Survival probability",
        main=paste(st,g,sprintf("log-rank p=%.3g",pp)))
      legend("topright",c("Low","High"),col=c("#4EAB8D","#C55286"),
        lty=1,bty="n",cex=.8)
      kmout[[length(kmout)+1L]]<-data.frame(cohort=st,gene=g,n=nrow(zz),
        median_cut=med,logrank_p=pp)
    }
    par(oldpar);dev.off()
    if(length(kmout))methods_csv(do.call(rbind,kmout),out,"08_Exploratory_KM_Logrank.csv")
  }
  methods_note(out,"NEW Cox implementation, not claimed to reproduce original 30 candidate selection.",
    "Input modules, cohort expression matrix and sample eligibility MUST be verified.",
    "If original rank, clinical adjustment, cohort-specific KM rule or BH family differs,",
    "published results must follow ORIGINAL code and verified input.")
  invisible(out)
}

# ===== M13 | Endpoint reporting / figures / software versions ===============
methods_M13 <- function() {
  out<-methods_dir("M13_Publication_Audit")
  keys<-c("ATLAS_RDS","ATLAS_META_CSV","MAC_RDS","AUC_RDS",
          "REGULON_RDS","SIX_STATE_RDS","SLINGSHOT_OLD",
          "PSEUDOBULK_MAP_CSV","PROGRAM_GENES_CSV","TF_MODULE_CSV")
  paths<-vapply(keys,function(k)PLOS[[k]],character(1))
  tab<-data.frame(asset=keys,path=paths,exists=file.exists(paths),
    bytes=ifelse(file.exists(paths),file.info(paths)$size,NA_real_))
  methods_csv(tab,out,"Methods_Inputs_Existence.csv")
  methods_note(out,"The ultimate E3 release must also include original data acquisition",
    "and processing scripts, package lockfile/sessionInfo, input checksums,",
    "true all-figure output maps, and legally shareable processed data.")
  invisible(tab)
}

# =============================================================================
# E3 REVIEW RECONSTRUCTION (2026-10-25) -- upgrades for M04/M04C/M06/M09/M11/M12
# These methods override the provisional implementations defined above.
# VERIFIED != IMPLEMENTED: new computations require comparison with frozen inputs.
# No method below overwrites pre-existing research objects or final figure files.
# =============================================================================

PLOS$FIG3_INPUT_DIR <- file.path(PLOS$ROOT,"E3_inputs","Fig3")
PLOS$PSEUDOBULK_MAP_CSV <- file.path(PLOS$FIG3_INPUT_DIR,
                                     "Fig3AB_verified_analysis_unit_map.csv")
PLOS$PROGRAM_GENES_CSV <- file.path(PLOS$FIG3_INPUT_DIR,
                                     "Fig3C_verified_program_gene_sets.csv")
PLOS$FIG5_INPUT_DIR <- file.path(PLOS$ROOT,"E3_inputs","Fig5")
PLOS$FIG5_MODULE_CSV <- file.path(PLOS$FIG5_INPUT_DIR,
                                  "Fig5_verified_three_module_membership.csv")
PLOS$FIG5_PARAMETERS <- file.path(PLOS$FIG5_INPUT_DIR,
                                  "Fig5_verified_network_parameters.csv")
PLOS$FIG6_FROZEN30 <- file.path(PLOS$SURV_DIR,"Fig6_verified_frozen_30_candidates.csv")

methods_write_template <- function(data, filename) {
  dir.create(dirname(filename),recursive=TRUE,showWarnings=FALSE)
  if (!file.exists(filename)) {
    utils::write.csv(data,filename,row.names=FALSE,na="")
    message("Created review template (NOT confirmed input): ",filename)
  }
  invisible(filename)
}
methods_grouped_counts <- function(ct,cells,groups,group_levels=NULL) {
  methods_need("Matrix")
  if(anyDuplicated(cells)||!all(cells %in% colnames(ct)))
    stop("Cell names must be unique and present in the count matrix")
  if(length(cells)!=length(groups)||anyNA(groups)||any(!nzchar(groups)))
    stop("Count aggregation: missing/invalid group identifiers")
  if(is.null(group_levels))group_levels<-unique(groups)
  g<-factor(groups,levels=group_levels)
  if(anyNA(g))stop("Count aggregation: unknown group label")
  if(any(table(g)==0L))stop("Count aggregation: unused/empty group")
  mat<-ct[,cells,drop=FALSE]
  design<-Matrix::sparse.model.matrix(~0+g)
  colnames(design)<-levels(g)
  if(nrow(design)!=ncol(mat))stop("Cell-by-group matrix alignment failure")
  ag<-mat%*%design
  colnames(ag)<-group_levels
  ag
}
methods_find_saved <- function(folders,regex,limit=300L) {
  folders<-folders[dir.exists(folders)]
  if(!length(folders))return(character())
  z<-unlist(lapply(folders,function(p)list.files(p,recursive=TRUE,
                         full.names=TRUE,pattern=regex,ignore.case=TRUE)),use.names=FALSE)
  sort(unique(head(z,limit)))
}
methods_assert_stage_samples <- function(meta,expect,source="") {
  m<-unique(meta[,c("sample_uid","stage"),drop=FALSE])
  if(anyDuplicated(m$sample_uid))stop(source,": same tissue ID assigned to two stages")
  actual<-table(factor(m$stage,levels=names(expect)))
  if(!identical(as.integer(actual),as.integer(expect)))
    stop(source,": expected ",paste(expect,collapse="/")," samples, found ",
         paste(actual,collapse="/"))
  invisible(actual)
}

# PREP A: explicit proposed 3/3/7/7 map for verification, NOT assumed donors.
methods_prepare_M04_mapping <- function() {
  out<-methods_dir("M04_Fig3AB_Pseudobulk_GSEA")
  meta<-methods_meta()
  if(!("cell" %in% names(meta)) || anyNA(meta$cell) || anyDuplicated(meta$cell))
    stop("Cell barcodes absent or duplicated in atlas metadata")
  req<-c("sample_uid","stage")
  if(!all(req %in% names(meta)))stop("Cannot construct mapping: missing stage/sample_uid")
  m<-meta[as.character(meta$stage) %in% PLOS$STAGES,
          c("cell","sample_uid","stage"),drop=FALSE]
  u<-unique(m[,c("sample_uid","stage")]);if(anyDuplicated(u$sample_uid))
    stop("Repeated sample UID across stages")
  sample_text<-sub("^.*::","",as.character(m$sample_uid))
  # The following GC collapse is an exploratory reconstruction of historic 7 GC
  # profiles from 9 samples and is NOT biological donor verification.
  base<-ifelse(as.character(m$stage)=="GC",
               sub("_t[0-9]+$","",sample_text),as.character(m$sample_uid))
  m$analysis_unit<-paste(m$stage,base,sep="::")
  z<-m[,c("cell","analysis_unit","stage","sample_uid")]
  summary<-unique(z[,c("analysis_unit","stage","sample_uid")])
  checks<-unique(summary[,c("analysis_unit","stage")])
  methods_csv(summary,out,"DRAFT_Fig3AB_Unit_By_Sample_MUST_VERIFY.csv")
  methods_csv(z,out,"DRAFT_Fig3AB_Cell_To_Unit_MUST_VERIFY.csv")
  methods_csv(data.frame(stage=PLOS$STAGES,n_profiles=as.integer(table(
    factor(checks$stage,levels=PLOS$STAGES)))),out,"DRAFT_Unit_Count_Audit.csv")
  methods_note(out,"DRAFT automatically inferred from sample labels, not independently verified",
      "Never claim GC prefixes encode verified donors without an original manifest.",
      "Verify sample-level pairing, then save a separate approved map at:",
      PLOS$PSEUDOBULK_MAP_CSV,
      "Approved file columns must include cell, analysis_unit, stage.")
  invisible(summary)
}
methods_M04 <- function() {
  out<-methods_dir("M04_Fig3AB_Pseudobulk_GSEA")
  if(!file.exists(PLOS$PSEUDOBULK_MAP_CSV)) {
    methods_prepare_M04_mapping()
    stop("M04 PAUSED: verify proposed sample/analysis-unit assignments first; ",
         "save independently approved mapping to ",PLOS$PSEUDOBULK_MAP_CSV,
         call.=FALSE)
  }
  methods_need(c("SeuratObject","Matrix","edgeR","limma","fgsea","msigdbr"))
  mp<-methods_read(PLOS$PSEUDOBULK_MAP_CSV,
                   c("cell","analysis_unit","stage"))
  if(anyNA(mp[,c("cell","analysis_unit","stage")])||anyDuplicated(mp$cell))
    stop("M04 mapping missing/duplicate cell information")
  if(!all(mp$stage %in% PLOS$STAGES))stop("Invalid stage in M04 mapping")
  u<-unique(mp[,c("analysis_unit","stage")])
  if(anyDuplicated(u$analysis_unit))stop("Analysis unit spans multiple stages")
  if(!identical(as.integer(table(factor(u$stage,levels=PLOS$STAGES))),
                c(3L,3L,7L,7L)))stop("Unexpected Fig3AB analysis units")
  x<-methods_load_seu();ct<-methods_assay_counts(x)
  if(!all(mp$cell %in% colnames(ct)))stop("Unknown cell in verified M04 map")
  # Prefer exact atlas cell ordering; avoids dense 81,779 x 20 transpose.
  mp<-mp[match(intersect(colnames(ct),mp$cell),mp$cell),,drop=FALSE]
  mat<-methods_grouped_counts(ct,mp$cell,mp$analysis_unit,u$analysis_unit)
  stages<-u$stage[match(colnames(mat),u$analysis_unit)]
  design<-stats::model.matrix(~0+factor(stages,levels=PLOS$STAGES))
  colnames(design)<-PLOS$STAGES
  y<-edgeR::DGEList(counts=mat,group=stages)
  keep<-edgeR::filterByExpr(y,design=design)
  y<-y[keep,,keep.lib.sizes=FALSE]
  y<-edgeR::calcNormFactors(y,method="TMM")
  vv<-limma::voom(y,design,plot=FALSE)
  fit<-limma::lmFit(vv,design)
  cc<-limma::makeContrasts(
    NAG=NAG-(CAG+IM+GC)/3,CAG=CAG-(NAG+IM+GC)/3,
    IM=IM-(NAG+CAG+GC)/3,GC=GC-(NAG+CAG+IM)/3,levels=design)
  fit<-limma::eBayes(limma::contrasts.fit(fit,cc))
  ms<-tryCatch(msigdbr::msigdbr(species="Homo sapiens",collection="H"),
               error=function(e)msigdbr::msigdbr(species="Homo sapiens",category="H"))
  hallmark<-lapply(split(ms$gene_symbol,ms$gs_name),unique)
  methods_csv(u,out,"01_Approved_Pseudobulk_Units.csv")
  methods_csv(data.frame(stage=PLOS$STAGES,profiles=as.integer(table(
    factor(u$stage,levels=PLOS$STAGES)))),out,"02_Stage_Unit_Counts.csv")
  for(st in PLOS$STAGES) {
    z<-limma::topTable(fit,coef=st,number=Inf,sort.by="none")
    z$gene<-rownames(z)
    methods_csv(z,out,paste0("Fig3A_DE_",st,".csv"))
    rank<-z$t;names(rank)<-z$gene
    rank<-rank[is.finite(rank)&!is.na(names(rank))&nzchar(names(rank))]
    rank<-sort(rank[!duplicated(names(rank))],decreasing=TRUE)
    fg<-as.data.frame(fgsea::fgseaMultilevel(pathways=hallmark,
           stats=rank,minSize=10,maxSize=500,eps=0))
    if("leadingEdge" %in% names(fg))
      fg$leadingEdge<-vapply(fg$leadingEdge,paste,character(1),collapse=";")
    methods_csv(fg,out,paste0("Fig3B_Hallmark_",st,".csv"))
  }
  methods_note(out,"NEW rerun based on independently approved 3/3/7/7 profile mapping.",
     "The profile IDs are not by definition independent patients.",
     "Interpret limma/voom P values as exploratory if biological units overlap.",
     paste("MSigDB/msigdbr package:",as.character(utils::packageVersion("msigdbr"))),
     "Historical gene-set release, exact filters and contrasts require comparison.")
  invisible(out)
}

# PREP B: program genes must be source-derived rather than guessed from labels.
methods_prepare_M04C <- function() {
  d<-methods_dir("M04C_Fig3C_Sample_Programs")
  methods_write_template(data.frame(program=c("Energy metabolism",
     "Inflammatory response","Intestinal differentiation"),gene=NA_character_,
     provenance="REPLACE with original program gene symbols"),
     PLOS$PROGRAM_GENES_CSV)
  saved<-methods_find_saved(c(file.path(PLOS$ATLAS_DIR,"Fig3_ABC_Revised"),
    PLOS$ATLAS_DIR),"(program|score|gene.set|Fig3C).*(csv|tsv|rds|R)$")
  methods_csv(data.frame(candidate_file=saved),d,"Candidate_Historical_Fig3C_Files.csv")
  methods_note(d,"Required: original exact gene set membership for three programs.",
    "File auto-created at: ",PLOS$PROGRAM_GENES_CSV,
    "DO NOT leave placeholder NA genes in the approved file.")
  invisible(d)
}
methods_M04C <- function() {
  out<-methods_dir("M04C_Fig3C_Sample_Programs")
  if(!file.exists(PLOS$PROGRAM_GENES_CSV)) {
    methods_prepare_M04C()
    stop("M04C PAUSED: supply the three original program gene sets",call.=FALSE)
  }
  gs<-methods_read(PLOS$PROGRAM_GENES_CSV,c("program","gene"))
  gs$gene<-trimws(as.character(gs$gene))
  if(anyNA(gs$gene)||any(!nzchar(gs$gene)))
    stop("M04C gene list contains placeholder/missing gene; verify first")
  if(length(unique(gs$program))!=3L)stop("Expected exact three program names")
  methods_need(c("SeuratObject","Matrix","edgeR","ggplot2"))
  x<-methods_load_seu();m<-x@meta.data
  if(!all(c("sample_uid","stage") %in% names(m)))stop("Atlas lacks sample/stage")
  keep<-which(as.character(m$stage) %in% PLOS$STAGES)
  sp<-data.frame(sample_uid=as.character(m$sample_uid[keep]),
                 stage=as.character(m$stage[keep]))
  methods_assert_stage_samples(sp,c(NAG=3L,CAG=3L,IM=7L,GC=9L),"Fig3C")
  uu<-unique(sp)
  mat<-methods_grouped_counts(methods_assay_counts(x),colnames(x)[keep],
                              sp$sample_uid,uu$sample_uid)
  y<-edgeR::DGEList(counts=mat)
  y<-edgeR::calcNormFactors(y,method="TMM")
  expr<-edgeR::cpm(y,log=TRUE,prior.count=1)
  genes<-unique(gs$gene)
  coverage<-data.frame(program=gs$program,gene=gs$gene,
                       detected=gs$gene %in% rownames(expr))
  methods_csv(coverage,out,"Fig3C_Gene_Coverage.csv")
  matched<-intersect(genes,rownames(expr))
  if(!length(matched))stop("No program genes represented")
  z<-t(scale(t(expr[matched,,drop=FALSE])))
  bad<-rownames(z)[rowSums(!is.finite(z))>0L]
  z<-z[rowSums(!is.finite(z))==0L,,drop=FALSE]
  if(!nrow(z))stop("No variable genes to score")
  by_program<-split(gs$gene,gs$program)
  res<-do.call(rbind,lapply(names(by_program),function(program_name){
    g<-by_program[[program_name]]
    k<-intersect(unique(g),rownames(z))
    if(!length(k))stop("Program has no valid variable genes")
    data.frame(sample_uid=colnames(z),program=program_name,
               score=colMeans(z[k,,drop=FALSE]),n_program_genes=length(k))
  }))
  res$stage<-uu$stage[match(res$sample_uid,uu$sample_uid)]
  methods_csv(data.frame(excluded_zero_variance=bad),out,"ZeroVariance_Genes.csv")
  methods_csv(res,out,"Fig3C_Tissue_Sample_Program_Scores.csv")
  gg<-ggplot2::ggplot(res,ggplot2::aes(stage,score))+
     ggplot2::geom_boxplot(outlier.shape=NA)+
     ggplot2::geom_point(position=ggplot2::position_jitter(width=.06,height=0))+
     ggplot2::facet_wrap(~program,scales="free_y")+
     ggplot2::theme_classic()+ggplot2::labs(y="Mean within-gene Z score",x=NULL)
  ggplot2::ggsave(file.path(out,"Fig3C_Score_Reconstruction.pdf"),gg,width=10,height=4)
  methods_note(out,"NEW rerun from frozen exact program gene list.",
    "9 GC tissue samples are NOT automatically 9 independent patients.",
    "Compare panel values and gene-wise standardization with the published frozen Fig3C.")
  invisible(res)
}

# M06: exact 3,101-candidate reprocessing, cluster comparison, 1,529 cell crosswalk.
methods_M06 <- function() {
  methods_need(c("Seurat","SeuratObject","ggplot2"))
  out<-methods_dir("M06_Macrophage_Seurat_Rerun")
  original<-methods_load_seu(PLOS$MAC_RDS)
  if(ncol(original)!=3101L)stop("M06 requires 3,101 original candidates")
  final<-methods_load_seu(PLOS$SIX_STATE_RDS)
  if(ncol(final)!=1529L)stop("Final six-state object must contain 1,529 cells")
  oldmeta<-original@meta.data
  origcl<-intersect(c("cluster_r1_0","seurat_clusters"),names(oldmeta))[1]
  if(is.na(origcl))stop("No frozen 17-cluster annotation available")
  oldlabel<-as.character(oldmeta[[origcl]])
  stcol<-intersect(c("fig4_state_display","Fig4_six_state"),
                   names(final@meta.data))[1]
  if(is.na(stcol))stop("No frozen six-state label")
  overlap<-intersect(colnames(original),colnames(final))
  if(!length(overlap))stop("Final and original cell barcodes do not overlap")
  if(length(overlap)!=ncol(final))stop("Not all final six-state cells trace to initial 3101")
  cross<-data.frame(cell=overlap,
     original_cluster=oldlabel[match(overlap,colnames(original))],
     final_state=as.character(final@meta.data[overlap,stcol]),
     stringsAsFactors=FALSE)
  methods_csv(cross,out,"01_Original_Cluster_to_Final_SixState_Crosswalk.csv")
  methods_csv(as.data.frame(table(original_cluster=cross$original_cluster,
                                        final_state=cross$final_state)),out,
              "02_Crosswalk_Counts.csv")
  methods_csv(as.data.frame(table(stage=as.character(final$stage),
     state=as.character(final@meta.data[[stcol]]))),out,
     "03_Final_SixState_by_Stage_Counts.csv")
  # Reprocessing deliberately starts from the frozen 3,101-cell input;
  # the original clustering labels are preserved for independent comparison.
  x<-original
  Seurat::DefaultAssay(x)<-"RNA"
  set.seed(PLOS$SEED)
  x<-Seurat::NormalizeData(x,verbose=FALSE)
  x<-Seurat::FindVariableFeatures(x,nfeatures=2500,verbose=FALSE)
  vf<-Seurat::VariableFeatures(x)
  x<-SeuratObject::`VariableFeatures<-`(x,value=vf[!grepl("^(MT-|RPL|RPS)",vf)])
  if(length(Seurat::VariableFeatures(x))<1000L)
    stop("Insufficient HVGs after exclusion of mitochondrial/ribosomal genes")
  x<-Seurat::ScaleData(x,verbose=FALSE)
  x<-Seurat::RunPCA(x,npcs=30,verbose=FALSE)
  x<-Seurat::FindNeighbors(x,dims=1:20,verbose=FALSE)
  x<-Seurat::FindClusters(x,resolution=1.0,verbose=FALSE)
  x<-Seurat::RunUMAP(x,dims=1:20,seed.use=PLOS$SEED,verbose=FALSE)
  z<-Seurat::FindAllMarkers(x,only.pos=FALSE,min.pct=.20,
                             logfc.threshold=.25)
  methods_csv(z,out,"04_Rerun_Cluster_Markers.csv")
  rt<-data.frame(cell=colnames(x),old=oldlabel[match(colnames(x),
                 colnames(original))],new=as.character(x$seurat_clusters))
  methods_csv(rt,out,"05_Rerun_vs_Frozen_Cluster_Assignments.csv")
  # ARI by contingency cells, label-invariant; no dependency on mclust.
  ari<-function(a,b) {
    t<-table(a,b);comb2<-function(y)sum(as.numeric(y)*(as.numeric(y)-1)/2)
    N<-sum(t);total<-N*(N-1)/2
    if(total==0)return(NA_real_)
    ab<-comb2(t);ax<-comb2(rowSums(t));bx<-comb2(colSums(t))
    expected<-ax*bx/total;mx<-(ax+bx)/2
    if(mx==expected)return(NA_real_)
    (ab-expected)/(mx-expected)
  }
  score<-data.frame(n_cells=ncol(x),n_original_clusters=length(unique(rt$old)),
    n_rerun_clusters=length(unique(rt$new)),adjusted_rand_index=ari(rt$old,rt$new),
    note="Independent reseeding; clusters need not have the same numeric IDs")
  methods_csv(score,out,"06_Rerun_ARI_vs_Frozen.csv")
  saveRDS(x,file.path(out,"Independent_HighRes_3101_Seurat.rds"))
  methods_note(out,
    "REBUILT 17-cluster candidate representation, compared with original by ARI.",
    "Final 1529-cell six-state labels were imported from the frozen annotated object.",
    "Cell types are not inferred from cluster indices; 150 NAT/exclusions must",
    "be traced with saved QC files in M07.",
    "This independent replay is NOT asserted identical to original 17 clusters.")
  print(score);invisible(out)
}

# M09: root sensitivity on the same frozen 1,529-cell embedding.
methods_M09 <- function() {
  methods_need(c("slingshot","SeuratObject"))
  out<-methods_dir("M09_Slingshot_Root_Sensitivity")
  x<-methods_load_seu(PLOS$SIX_STATE_RDS)
  if(ncol(x)!=1529L)stop("Expected final 1,529-cell six-state subset")
  state_col<-intersect(c("fig4_state_display","Fig4_six_state"),
                        names(x@meta.data))[1]
  if(is.na(state_col))stop("Missing verified six-state identities")
  state<-as.character(x@meta.data[[state_col]])
  if(anyNA(state))stop("Unannotated cells in frozen Slingshot subset")
  reduction<-intersect(c("umap","umap_highres","umap_final"),
                        names(x@reductions))[1]
  if(is.na(reduction))stop("Frozen 1,529-cell UMAP embedding absent")
  coords<-as.matrix(SeuratObject::Embeddings(x[[reduction]])[,1:2,drop=FALSE])
  if(!setequal(rownames(coords),colnames(x)))stop("Coordinate/cell mismatch")
  coords<-coords[colnames(x),,drop=FALSE]
  fc<-unique(state[grepl("FCN1",state)]);
  mixed<-unique(state[grepl("HLA-II|GAST",state)])
  if(length(fc)!=1L||length(mixed)!=1L)stop("Cannot resolve both roots")
  cases<-list(FCN1_six=list(root=fc,keep=rep(TRUE,length(state))),
    HLA_GAST_six=list(root=mixed,keep=rep(TRUE,length(state))),
    FCN1_five=list(root=fc,keep=state!=mixed))
  fits<-list();allpaths<-list();allpt<-list()
  for(nm in names(cases)) {
    conf<-cases[[nm]];idx<-which(conf$keep)
    f<-slingshot::slingshot(coords[idx,,drop=FALSE],
       clusterLabels=state[idx],start.clus=conf$root)
    fits[[nm]]<-f
    saveRDS(f,file.path(out,paste0("Slingshot_",nm,".rds")))
    paths<-slingshot::slingLineages(f)
    if(!length(paths))stop("Slingshot produced no paths in ",nm)
    allpaths[[nm]]<-data.frame(scenario=nm,root=conf$root,
       lineage=seq_along(paths),
       path=vapply(paths,paste,character(1),collapse=" -> "))
    pt<-as.matrix(slingshot::slingPseudotime(f,na=TRUE))
    colnames(pt)<-paste0("L",seq_len(ncol(pt)))
    allpt[[nm]]<-data.frame(cell=colnames(x)[idx],scenario=nm,pt,
                            check.names=FALSE)
  }
  methods_csv(do.call(rbind,allpaths),out,"01_Root_Lineage_Topologies.csv")
  methods_csv(do.call(rbind,allpt),out,"02_PerCell_Lineage_Pseudotime.csv")
  # Compare common six-state roots with one another by lineage sequence,
  # not by claiming cell development follows the same physical direction.
  p<-do.call(rbind,allpaths)
  tops<-lapply(split(p$path,p$scenario),unique)
  jj<-function(a,b){length(intersect(a,b))/length(union(a,b))}
  comp<-data.frame(comparison=c("FCN1_six_vs_HLA_GAST_six",
      "FCN1_six_vs_FCN1_five"),
    exact_lineage_sequence_Jaccard=c(jj(tops$FCN1_six,tops$HLA_GAST_six),
                                    jj(tops$FCN1_six,tops$FCN1_five)))
  methods_csv(comp,out,"03_Lineage_Topology_Sensitivity.csv")
  methods_csv(data.frame(stage=as.character(x$stage),state=state,
     cell=colnames(x)),out,"04_SixState_Stage_Cell_Map.csv")
  methods_note(out,"Root=FCN1 is an analytical hypothesis, not a demonstrated biological ancestor.",
    "Root=HLA-II/GAST is a sensitivity test; enrichment in early stages does not",
    "prove developmental priority for lineage-ambiguous cells.",
    paste("Frozen embedding:",reduction),
    "Fig4D reference requires separately checking original saved Slingshot fit.",
    "Exact ordered path Jaccard is a strict descriptor, not a formal test of stability.")
  print(comp);invisible(out)
}

# M11: matched final module provenance + reconstructed target overlap and GO.
methods_prepare_M11 <- function() {
  d<-methods_dir("M11_TF_Modules_GO")
  paths<-methods_find_saved(c(PLOS$TF_DIR,
      file.path(PLOS$FIG4_DIR,"Revision_Method_Audit")),
      "(network|module|regulon|overlap|GO|enrich|nodes|edges).*(csv|tsv|xlsx|rds|R)$",500L)
  methods_csv(data.frame(candidate_file=paths),d,
              "00_Locate_Saved_Network_Module_GO_Files.csv")
  methods_write_template(data.frame(TF=NA_character_,module=NA_character_,
    regulon=NA_character_,source_file="REPLACE with frozen final membership"),
    PLOS$FIG5_MODULE_CSV)
  methods_write_template(data.frame(parameter=c("network_scope",
    "edge_min_shared_targets","edge_min_jaccard","community_method",
    "module_selection_rule","regulon_family","GO_background_source"),
    verified_value=NA_character_),PLOS$FIG5_PARAMETERS)
  methods_note(d,"Original full vs top100 networks had different module partitions.",
    "The required FINAL selected-three-module membership cannot be inferred",
    "from numeric module labels alone; please verify with final Fig5 D/E analysis.")
  invisible(d)
}
methods_regulon_genes <- function(x) {
  if(is.character(x))return(unique(x[!is.na(x)&nzchar(x)]))
  if(is.numeric(x)&&!is.null(names(x)))return(unique(names(x)))
  if(is.list(x)&&!is.null(x$genes)) {
    if(is.numeric(x$genes)&&!is.null(names(x$genes)))
      return(unique(names(x$genes)))
    return(unique(as.character(x$genes)))
  }
  if(is.list(x)) {
    y<-unlist(x,use.names=FALSE)
    if(is.character(y))return(unique(y[!is.na(y)&nzchar(y)]))
  }
  character()
}
methods_M11 <- function() {
  out<-methods_dir("M11_TF_Modules_GO")
  reg<-readRDS(methods_file(PLOS$REGULON_RDS,"saved SCENIC regulons"))
  if(!is.list(reg)||is.null(names(reg)))stop("Expected named regulon list")
  methods_csv(data.frame(regulon=names(reg),
    target_count=vapply(reg,function(z)length(methods_regulon_genes(z)),integer(1))),
    out,"01_All_Saved_Regulons_Target_Count.csv")
  if(!file.exists(PLOS$FIG5_MODULE_CSV)) {
    methods_prepare_M11()
    stop("M11 PAUSED: final module-to-regulon assignments must be verified at ",
         PLOS$FIG5_MODULE_CSV,call.=FALSE)
  }
  m<-methods_read(PLOS$FIG5_MODULE_CSV,c("TF","module","regulon"))
  if(anyNA(m[,c("TF","module","regulon")])||anyDuplicated(m$TF))
    stop("Module input has placeholder/duplicate TFs")
  if(nrow(m)!=26L||length(unique(m$module))!=3L)
    stop("Expected exactly 26 TFs in exactly 3 frozen modules")
  if(!all(m$regulon %in% names(reg)))stop("Some frozen regulon identifiers absent")
  tf_parse<-sub("(_extended)? \\([0-9]+g\\)$","",m$regulon)
  if(!all(tf_parse==m$TF))stop("TF/regulon identifiers disagree")
  methods_csv(m,out,"02_Verified_26TF_Module_Identity.csv")
  genes<-setNames(lapply(m$regulon,function(n)methods_regulon_genes(reg[[n]])),m$TF)
  if(any(lengths(genes)==0L))stop("Frozen modules contain empty target list")
  gene_module<-do.call(rbind,lapply(seq_len(nrow(m)),function(i){
    data.frame(module=m$module[i],TF=m$TF[i],regulon=m$regulon[i],
               gene=genes[[i]],stringsAsFactors=FALSE)
  }))
  methods_csv(gene_module,out,"03_Verified_Module_TF_Target_Edges.csv")
  modsets<-lapply(split(m$TF,m$module),function(tf)
     unique(unlist(genes[tf],use.names=FALSE)))
  summary<-data.frame(module=names(modsets),
     n_TFs=as.integer(table(m$module)[names(modsets)]),
     n_unique_targets=lengths(modsets))
  methods_csv(summary,out,"04_Module_Target_Unions.csv")
  unique_genes<-unique(gene_module$gene)
  audit<-data.frame(n_TFs=nrow(m),n_modules=length(modsets),
    n_unique_predicted_targets=length(unique_genes),
    matches_manuscript_4387=length(unique_genes)==4387L)
  methods_csv(audit,out,"05_Manuscript_Target_Count_Check.csv")
  # Pairwise edges are exactly derived from the 26 frozen TF target sets.
  two<-utils::combn(m$TF,2L,simplify=FALSE)
  ed<-do.call(rbind,lapply(two,function(p){
    k<-length(intersect(genes[[p[1]]],genes[[p[2]]]))
    n<-length(union(genes[[p[1]]],genes[[p[2]]]))
    data.frame(from=p[1],to=p[2],shared_targets=k,
      jaccard=if(n)k/n else 0,module_from=m$module[match(p[1],m$TF)],
      module_to=m$module[match(p[2],m$TF)])
  }))
  methods_csv(ed,out,"06_Frozen26TF_AllPair_Target_Overlap.csv")
  # Publication GO comparison is meaningful only with the EXACT original
  # expressed-gene background and original module membership provided.
  go_bg<-file.path(PLOS$FIG5_INPUT_DIR,"Fig5_verified_GO_background_genes.csv")
  if(!file.exists(go_bg)) {
    methods_write_template(data.frame(gene=NA_character_),go_bg)
    methods_note(out,"GO analysis not executed: exact original gene background missing.",
     "Saved target/network overlap is reproducible from verified module assignment.",
     "A 26 TF / 4387 gene match is a numerical check, not proof that the original",
     "network community detection or final figure was regenerated.")
    message("M11 GO paused pending verified enrichment background")
    return(invisible(audit))
  }
  bg<-methods_read(go_bg,"gene")
  if(anyNA(bg$gene)||!length(bg$gene))stop("Unverified/empty GO universe")
  methods_need(c("clusterProfiler","org.Hs.eg.db"))
  keymap<-function(symbols) {
    z<-suppressMessages(clusterProfiler::bitr(unique(symbols),
       fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db::org.Hs.eg.db))
    unique(z$ENTREZID)
  }
  background<-keymap(bg$gene)
  for(mod in names(modsets)) {
    gi<-keymap(modsets[[mod]])
    if(length(gi)<10L)stop("Fewer than 10 target Entrez IDs for module ",mod)
    en<-clusterProfiler::enrichGO(gene=gi,universe=background,
      OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType="ENTREZID",ont="BP",
      pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,readable=TRUE)
    if(!is.null(en)&&nrow(as.data.frame(en)))
      methods_csv(as.data.frame(en),out,
                  paste0("07_Recomputed_GO_BP_Module_",gsub("[^A-Za-z0-9]","_",mod),".csv"))
  }
  methods_note(out,"GO results recomputed with a verified gene background.",
    "Check term selection/global-FDR derivation against the frozen Fig5E plot.",
    "Do not substitute per-module BH for a reported globally adjusted FDR.")
  invisible(audit)
}

# M12: freeze discovered candidates, audit exact 30->27->9->7 selection.
methods_prepare_M12 <- function() {
  d<-methods_dir("M12_TCGA_External_Cox")
  src<-methods_find_saved(c(PLOS$SURV_DIR,PLOS$ROOT),
       "(cox|survival|TCGA|62254|84437|external|clinical|candidate).*(csv|tsv|xlsx|rds|R)$",500L)
  methods_csv(data.frame(candidate_file=src),d,
     "00_Saved_Cox_Clinical_Expression_File_Inventory.csv")
  methods_write_template(data.frame(module=NA_character_,gene=NA_character_),
     PLOS$FIG6_FROZEN30)
  for(cohort in c("TCGA","GSE62254","GSE84437")) {
    methods_write_template(data.frame(sample_id=NA_character_,os_months=NA_real_,
      event=NA_integer_),file.path(PLOS$SURV_DIR,
                                  paste0(cohort,"_clinical.csv")))
    methods_write_template(data.frame(sample_id=NA_character_,
      IGFBP7=NA_real_),file.path(PLOS$SURV_DIR,
                               paste0(cohort,"_expression.csv")))
  }
  methods_note(d,"CSV layouts are templates, NOT sourced clinical cohorts.",
    "DO NOT use templates without replacing with verified original data.",
    "Frozen 30 candidate list must match final original Fig6A, not a fresh p-rank guess.",
    "Cohort time unit must be months and event coding must be confirmed.")
  invisible(d)
}
methods_M12 <- function() {
  methods_need("survival")
  out<-methods_dir("M12_TCGA_External_Cox")
  cohorts<-c("TCGA","GSE62254","GSE84437")
  need<-c(PLOS$FIG6_FROZEN30,unlist(lapply(cohorts,function(z)
     file.path(PLOS$SURV_DIR,paste0(z,c("_clinical.csv","_expression.csv"))))))
  if(!all(file.exists(need))) {
    methods_prepare_M12()
    stop("M12 PAUSED: clinical/expression cohorts and/or frozen 30 genes absent; ",
      "see templates and inventory under ",out,call.=FALSE)
  }
  fixed<-methods_read(PLOS$FIG6_FROZEN30,c("module","gene"))
  if(anyNA(fixed[,c("module","gene")])||nrow(fixed)!=30L||
     anyDuplicated(fixed$gene)||length(unique(fixed$module))!=3L||
     any(as.integer(table(fixed$module))!=10L))
    stop("Frozen selection must be exactly 10 unique genes from each of 3 modules")
  methods_csv(fixed,out,"01_Verified_Frozen_30_Candidates.csv")
  load_cohort<-function(id){
    cl<-methods_read(file.path(PLOS$SURV_DIR,paste0(id,"_clinical.csv")),
       c("sample_id","os_months","event"))
    ex<-methods_read(file.path(PLOS$SURV_DIR,paste0(id,"_expression.csv")),"sample_id")
    if(anyNA(cl$sample_id)||anyNA(ex$sample_id)||anyDuplicated(cl$sample_id)||
       anyDuplicated(ex$sample_id))stop("Missing or duplicate cohort sample IDs in ",id)
    if(!is.numeric(cl$os_months))stop("OS must be numeric MONTHS in ",id)
    if(!all(is.na(cl$event)|cl$event %in% c(0,1)))
      stop("Event must be coded 0/1 in ",id)
    joined<-merge(cl,ex,by="sample_id")
    good<-is.finite(joined$os_months)&joined$os_months>0&
          !is.na(joined$event)&joined$event %in% c(0,1)
    aud<-data.frame(cohort=id,clinical_rows=nrow(cl),expression_rows=nrow(ex),
      matched_samples=nrow(joined),invalid_OS_or_event=sum(!good),
      eligible_samples=sum(good),events=sum(joined$event[good]))
    list(data=joined[good,,drop=FALSE],audit=aud)
  }
  loads<-setNames(lapply(cohorts,load_cohort),cohorts)
  audit<-do.call(rbind,lapply(loads,`[[`,"audit"))
  methods_csv(audit,out,"02_Eligibility_And_Event_Audit.csv")
  expect<-c(TCGA=370L,GSE62254=300L,GSE84437=431L)
  event_expect<-c(TCGA=147L,GSE62254=152L,GSE84437=207L)
  if(!identical(setNames(as.integer(audit$eligible_samples),audit$cohort),expect)||
     !identical(setNames(as.integer(audit$events),audit$cohort),event_expect))
     stop("COHORT AUDIT FAILED: expected 370/300/431 with 147/152/207 events;",
       " inspect the input origin and inclusion criteria, don't force counts")
  dat<-lapply(loads,`[[`,"data")
  coverage<-expand.grid(gene=fixed$gene,cohort=cohorts,stringsAsFactors=FALSE)
  coverage$present<-mapply(function(g,id) g %in% names(dat[[id]])&&
      sum(is.finite(suppressWarnings(as.numeric(dat[[id]][[g]]))))>0L,
      coverage$gene,coverage$cohort)
  methods_csv(coverage,out,"03_Frozen_Gene_Coverage.csv")
  measurable<-vapply(cohorts,function(st)sum(coverage$present[coverage$cohort==st]),integer(1))
  shared<-sum(vapply(fixed$gene,function(g)all(coverage$present[coverage$gene==g &
                 coverage$cohort %in% c("GSE62254","GSE84437")]),logical(1)))
  cv<-data.frame(cohort=cohorts,measurable_genes=as.integer(measurable),
    shared_external_measurable=shared)
  methods_csv(cv,out,"04_Gene_Coverage_28_29_27_Check.csv")
  if(measurable["GSE62254"]!=28L||measurable["GSE84437"]!=29L||shared!=27L)
     stop("External panel coverage differs from manuscript 28/29/27; verify identifiers")
  cox_gene<-function(d,g,vars=character()) {
    if(!(g %in% names(d)))return(NULL)
    x<-suppressWarnings(as.numeric(d[[g]]))
    yy<-data.frame(time=d$os_months,event=d$event,expression=x)
    if(length(vars)) {
      if(!all(vars %in% names(d)))stop("Missing verified clinical covariates")
      if(any(make.names(vars)!=vars))stop("Use syntactically valid covariate column names")
      for(v in vars) yy[[v]]<-d[[v]]
    }
    yy<-yy[stats::complete.cases(yy),,drop=FALSE]
    if(nrow(yy)<10||sum(yy$event)<4||length(unique(yy$expression))<2)return(NULL)
    formula<-stats::as.formula(paste("survival::Surv(time,event) ~ expression",
      if(length(vars))paste("+",paste(vars,collapse=" + ")) else ""))
    fit<-tryCatch(survival::coxph(formula,data=yy),error=function(e)NULL)
    if(is.null(fit))return(NULL)
    ss<-summary(fit)
    if(!("expression" %in% rownames(ss$coefficients)))return(NULL)
    z<-ss$coefficients["expression",];ci<-ss$conf.int["expression",]
    data.frame(gene=g,n=nrow(yy),events=sum(yy$event),HR=unname(ci["exp(coef)"]),
       HR_lower95=unname(ci["lower .95"]),HR_upper95=unname(ci["upper .95"]),
       p_value=unname(z["Pr(>|z|)"]))
  }
  uni<-do.call(rbind,lapply(cohorts,function(id){
    z<-do.call(rbind,lapply(fixed$gene,function(g)cox_gene(dat[[id]],g)))
    if(is.null(z))stop("No Cox fits in cohort ",id)
    z$cohort<-id;z
  }))
  methods_csv(uni,out,"05_ThreeCohort_FrozenPanel_Univariable_Cox.csv")
  ext<-uni[uni$cohort!="TCGA",,drop=FALSE]
  for(id in c("GSE62254","GSE84437")) {
    idx<-which(ext$cohort==id)
    ext$BH_evaluable[idx]<-stats::p.adjust(ext$p_value[idx],method="BH")
    ext$BH_fixed30[idx]<-stats::p.adjust(ext$p_value[idx],method="BH",n=30L)
  }
  methods_csv(ext,out,"06_External_Cox_BH_Fixed30_vs_Evaluable.csv")
  result<-lapply(fixed$gene,function(g){
    a<-uni[uni$gene==g,,drop=FALSE]
    e<-ext[ext$gene==g,,drop=FALSE]
    direction<-nrow(a)==3L&&length(unique(sign(log(a$HR))))==1L
    nominal<-nrow(a)==3L&&all(a$p_value<0.05)
    data.frame(gene=g,module=fixed$module[match(g,fixed$gene)],
      same_HR_direction_all_three=direction,
      nominal_P_lt_0_05_all_three=nominal,
      external_FDR_lt_0_05_fixed30=nrow(e)==2L&&all(e$BH_fixed30<0.05),
      external_FDR_lt_0_05_evaluable=nrow(e)==2L&&all(e$BH_evaluable<0.05))
  })
  final<-do.call(rbind,result)
  final$nine_criterion<-final$same_HR_direction_all_three&
                         final$nominal_P_lt_0_05_all_three
  final$seven_criterion<-final$nine_criterion&final$external_FDR_lt_0_05_fixed30
  methods_csv(final,out,"07_CrossCohort_Consistency_30_to_9_to_7.csv")
  n9<-sum(final$nine_criterion);n7<-sum(final$seven_criterion)
  methods_csv(data.frame(frozen=30,shared_external=shared,
      n_crosscohort_nominal_direction=n9,n_external_fdr_fixed30=n7,
      matches_9_and_7=n9==9L&&n7==7L),out,"08_Manuscript_Validation_Check.csv")
  if(n9!=9L||n7!=7L)warning("Computed selection does not match manuscript 9/7;",
    " compare original Cox code, platform transformation and BH family")
  adj<-list()
  for(id in cohorts) {
    vars<-PLOS$COVARIATE_LIST[[id]]
    if(!length(vars))next
    z<-do.call(rbind,lapply(fixed$gene,function(g)cox_gene(dat[[id]],g,vars)))
    if(!is.null(z)){z$cohort<-id;z$adjustment<-paste(vars,collapse=";");adj[[id]]<-z}
  }
  if(length(adj))methods_csv(do.call(rbind,adj),out,
                      "09_Verified_Covariate_Adjusted_Cox.csv")
  # KM logrank: explicit cohort-specific median, separate inference from Cox.
  reps<-c("IGFBP7","SPARC","IKBKE")
  km<-list()
  for(id in cohorts)for(g in reps){
    d<-dat[[id]];if(!g %in% names(d))next
    x<-suppressWarnings(as.numeric(d[[g]]))
    valid<-is.finite(x)
    if(sum(valid)<10L)next
    y<-data.frame(time=d$os_months[valid],event=d$event[valid],x=x[valid]);
    cut<-stats::median(y$x)
    y$group<-factor(ifelse(y$x>cut,"High","Low"))
    if(length(unique(y$group))!=2L)next
    zz<-survival::survdiff(survival::Surv(time,event)~group,data=y)
    km[[paste(id,g)]]<-data.frame(cohort=id,gene=g,n=nrow(y),
      median_expression=cut,logrank_p=stats::pchisq(zz$chisq,df=1,
                                                  lower.tail=FALSE))
  }
  if(length(km))methods_csv(do.call(rbind,km),out,"10_Representative_KM_Logrank.csv")
  methods_note(out,"Reanalysis uses independently frozen 30-gene panel from original Fig6.",
   "Eligibility, events, gene coverage, univariable Cox and BH scopes audited.",
   "No clinical adjustment unless verified covariate names are supplied.",
   "No exact discovery 4387-gene reranking claimed until original rank source is verified.",
   "HR is per one unit of supplied expression; platform scaling must match original.",
   "Genes passing Cox and genes passing median-split logrank need not coincide.")
  print(c(n9=n9,n7=n7));invisible(out)
}

# Audit helpers: call without expensive analysis; no information fabricated.
methods_missing_inputs <- function() {
  names<-c("Fig3AB profile mapping","Fig3C exact gene sets",
           "Frozen 3101 candidate Seurat","Frozen 1529 six-state Seurat",
           "Fig5 final three module IDs","Fig5 final network parameters",
           "Fig6 original frozen 30-gene list",
           "TCGA clinical","TCGA expression",
           "GSE62254 clinical","GSE62254 expression",
           "GSE84437 clinical","GSE84437 expression")
  paths<-c(PLOS$PSEUDOBULK_MAP_CSV,PLOS$PROGRAM_GENES_CSV,
           PLOS$MAC_RDS,PLOS$SIX_STATE_RDS,PLOS$FIG5_MODULE_CSV,
           PLOS$FIG5_PARAMETERS,PLOS$FIG6_FROZEN30,
           unlist(lapply(c("TCGA","GSE62254","GSE84437"),function(id)
             file.path(PLOS$SURV_DIR,paste0(id,c("_clinical.csv","_expression.csv"))))))
  z<-data.frame(required_input=names,file_path=paths,exists=file.exists(paths))
  print(z,row.names=FALSE)
  invisible(z)
}

METHODS <- list(
 M00=list(title="Software, execution environment",status="NEW",fun=methods_capture_env),
 M01=list(title="Data sources / original atlas metadata",status="SAVED",fun=methods_M01),
 M02=list(title="Seurat/Harmony annotated atlas inspection",status="SAVED",fun=methods_M02),
 M03=list(title="Fig1 twelve cell populations and fractions",status="NEW",fun=methods_M03),
 M04=list(title="Fig3AB pseudobulk + Hallmark GSEA",status="CODE_READY_INPUT_UNVERIFIED",fun=methods_M04),
 M04C=list(title="Fig3C 3/3/7/9 sample program scores",status="CODE_READY_INPUT_UNVERIFIED",fun=methods_M04C),
 M05=list(title="Fig2 CellChat saved networks and centrality",status="SAVED",fun=methods_M05),
 M06=list(title="Fig4 high-resolution Seurat reclustering",status="CODE_READY_NEEDS_COMPARISON",fun=methods_M06),
 M07=list(title="Fig4 scDblFinder DecontX saved QC",status="SAVED",fun=methods_M07),
 M08=list(title="Fig3 bulk CIBERSORT saved-results audit",status="SAVED",fun=methods_M08),
 M09=list(title="Fig4 Slingshot two-root sensitivity",status="CODE_READY_NEEDS_COMPARISON",fun=methods_M09),
 M10=list(title="Fig5 SCENIC 326 x 1529 AUCell",status="SAVED+NEW",fun=methods_M10),
 M11=list(title="Fig5 co-target network / GO provenance",status="CODE_READY_INPUT_UNVERIFIED",fun=methods_M11),
 M12=list(title="Fig6 TCGA Cox + two external validations",status="CODE_READY_INPUT_UNVERIFIED",fun=methods_M12),
 M13=list(title="Final source and version audit",status="NEW",fun=methods_M13)
)
methods_status <- function(){
  print(data.frame(module=names(METHODS),
    content=vapply(METHODS,`[[`,character(1),"title"),
    status=vapply(METHODS,`[[`,character(1),"status")),row.names=FALSE)
  cat("\nCODE_READY != FIGURE_REPRODUCED. Run methods_missing_inputs() first.\n")
  cat("Run one at a time: methods_run('M01')\n")
  cat("Outputs only in new directory: ",PLOS$OUT,"\n",sep="")
  invisible(METHODS)
}
methods_run <- function(which){
  stopifnot(is.character(which),length(which)>=1L)
  unknown<-setdiff(which,names(METHODS))
  if(length(unknown))stop("Unknown module(s): ",paste(unknown,collapse=", "))
  set.seed(PLOS$SEED)
  for(nm in which){
    cat("\n========== START ",nm," : ",METHODS[[nm]]$title," ==========\n",sep="")
    t0<-Sys.time()
    # Fail loudly and stop. No silent module skipping or fabricated outputs.
    result<-METHODS[[nm]]$fun()
    cat("========== DONE ",nm," in ",round(as.numeric(difftime(Sys.time(),t0,units="secs")),1)," sec ==========\n",sep="")
  }
  invisible(PLOS$OUT)
}

# =============================================================================
# ARCHIVED ACTUAL PROJECT R SCRIPTS, embedded below as OPTIONAL FUNCTIONS.
# These are exact recovered historical versions with hard-coded original paths.
# They do NOT execute on source(). Their output directories may be old and
# may overwrite old reports! Do NOT invoke unless paths are manually reviewed.
# These are NOT automatically attributed to current published final figures.
# =============================================================================

# ---- legacy_01_GeneUniverse_GateA : Gene audit preliminary ----
legacy_01_GeneUniverse_GateA <- function() {
# ===================================================================
# GSE150290 Gate A / Phase 3: Gene-universe and marker-coverage audit
# Standalone script | Base R only | Windows R 4.3.x
# No Seurat merging, normalization, clustering, or modifications to a1
# ===================================================================

OUT_DIR  <- "E:/download/GSE150290_audit"
DATA_DIR <- file.path(OUT_DIR, "processed_preview")
QC_FILE  <- file.path(OUT_DIR, "05_all52_sample_QC_summary.csv")

stopifnot(dir.exists(DATA_DIR), file.exists(QC_FILE))
qc <- read.csv(QC_FILE, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(qc) == 52L, !anyDuplicated(qc$sample))

files <- list.files(DATA_DIR, pattern = "\\.txt\\.gz$", recursive = TRUE,
                    full.names = TRUE)
stopifnot(length(files) == 52L)
fn <- basename(files)
gsm <- sub("^(GSM[0-9]+)_.*$", "\\1", fn)
sample <- sub("^GSM[0-9]+_(Pat[0-9]+-[A-Za-z]).*$", "\\1", fn)
if (anyDuplicated(sample) || !setequal(sample, qc$sample) ||
    any(!grepl("^GSM[0-9]+$", gsm)) ||
    any(!grepl("^Pat[0-9]+-[A-Za-z]$", sample))) {
  stop("Processed files do not match the 52-sample QC manifest")
}
idx <- match(qc$sample, sample)
files <- files[idx]

# Presence is not cell-identity proof. These markers are only a feasibility screen.
markers <- list(
  Macrophage = c("LST1", "TYROBP", "FCER1G", "AIF1", "CD68", "C1QA",
                 "C1QB", "MS4A7", "CSF1R", "CD163"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  T_cell = c("CD3D", "CD3E", "TRAC"),
  B_cell = c("MS4A1", "CD79A"),
  Endothelial = c("PECAM1", "VWF", "KDR"),
  Fibroblast = c("COL1A1", "DCN", "LUM"),
  Key_CellChat = c("MIF", "CD74", "CD44", "CXCR4", "TNFSF13B", "TNFRSF17"),
  Key_Survival = c("BASP1", "MAP4K4", "NR1D1", "SERPINE1", "KBTBD2", "VAT1")
)
marker_keys <- unique(unlist(markers, use.names = FALSE))
marker_category <- vapply(marker_keys, function(g) {
  paste(names(markers)[vapply(markers, function(m) g %in% m, logical(1))],
        collapse = ";")
}, character(1))

scan_one <- function(path, expected_cells, wanted, chunk_size = 250L) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  head <- readLines(con, n = 1L, warn = FALSE)
  if (length(head) != 1L) stop("Unreadable header: ", path)
  n_cells <- length(strsplit(head, "\t", fixed = TRUE)[[1]])
  if (n_cells != expected_cells) stop("Cell columns mismatch: ", basename(path))
  genes <- character()
  observed <- list()
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    g <- sub("\t.*$", "", lines)
    genes <- c(genes, g)
    positions <- which(g %in% wanted)
    if (length(positions)) {
      for (i in positions) {
        parts <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
        if (length(parts) != n_cells + 1L) {
          stop("Unexpected row width for ", g[i], " in ", basename(path))
        }
        z <- suppressWarnings(as.numeric(parts[-1L]))
        if (any(!is.finite(z))) stop("Non-finite counts for ", g[i])
        if (!is.null(observed[[g[i]]])) stop("Duplicate marker gene: ", g[i])
        observed[[g[i]]] <- c(n_cells_expr = sum(z > 0),
                             detection_pct = 100 * mean(z > 0),
                             mean_count = mean(z))
      }
    }
  }
  list(genes = genes, markers = observed, n_cells = n_cells)
}

all_gene_sets <- vector("list", nrow(qc))
marker_rows <- vector("list", nrow(qc))
for (i in seq_len(nrow(qc))) {
  cat(sprintf("[%02d/52] %-8s scanning gene universe and selected markers...\n",
              i, qc$sample[i]))
  result <- scan_one(files[i], qc$n_cells_matrix[i], marker_keys)
  all_gene_sets[[i]] <- result$genes
  if (anyDuplicated(result$genes)) stop("Duplicated gene names: ", qc$sample[i])
  if (length(result$genes) != qc$n_genes[i]) {
    stop("Reported gene count mismatch: ", qc$sample[i])
  }
  marker_rows[[i]] <- do.call(rbind, lapply(marker_keys, function(g) {
    z <- result$markers[[g]]
    data.frame(sample = qc$sample[i], GSM = qc$GSM[i],
               gene = g, category = marker_category[[g]],
               detected_gene = !is.null(z),
               n_cells_expr = if (is.null(z)) NA_integer_ else unname(z["n_cells_expr"]),
               detection_pct = if (is.null(z)) NA_real_ else unname(z["detection_pct"]),
               mean_count = if (is.null(z)) NA_real_ else unname(z["mean_count"]),
               stringsAsFactors = FALSE)
  }))
  rm(result)
}
names(all_gene_sets) <- qc$sample

ref_early <- all_gene_sets[["Pat25-A"]]
ref_other <- all_gene_sets[["Pat01-A"]]
shared_all <- Reduce(intersect, all_gene_sets)
union_all <- Reduce(union, all_gene_sets)

coverage <- data.frame(
  sample = qc$sample,
  patient = qc$patient,
  candidate_group = qc$group_candidate,
  n_genes = lengths(all_gene_sets),
  overlap_with_Pat01 = vapply(all_gene_sets, function(x) length(intersect(x, ref_other)), integer(1)),
  overlap_with_Pat25 = vapply(all_gene_sets, function(x) length(intersect(x, ref_early)), integer(1)),
  identical_set_to_Pat01 = vapply(all_gene_sets, function(x) setequal(x, ref_other), logical(1)),
  identical_set_to_Pat25 = vapply(all_gene_sets, function(x) setequal(x, ref_early), logical(1)),
  stringsAsFactors = FALSE
)
marker_table <- do.call(rbind, marker_rows)
rownames(marker_table) <- NULL

# A useful first-pass per-sample panel completeness report; not a cell annotation.
panel_summary <- do.call(rbind, lapply(names(markers), function(category) {
  gs <- markers[[category]]
  d <- subset(marker_table, gene %in% gs)
  x <- aggregate(detected_gene ~ sample, d, sum)
  x$category <- category
  x$n_panel_genes <- length(gs)
  names(x)[2] <- "n_detected_genes"
  x
}))
panel_summary <- merge(panel_summary, qc[, c("sample", "patient", "group_candidate", "n_cells_pass_qc")],
                       by = "sample", all.x = TRUE, sort = FALSE)

# Save all outputs in a separate phase 3 folder; do not overwrite phases 1-2.
P3 <- file.path(OUT_DIR, "GateA_Phase3_gene_audit")
dir.create(P3, recursive = TRUE, showWarnings = FALSE)
write.csv(coverage, file.path(P3, "09_gene_universe_by_sample.csv"), row.names = FALSE)
write.csv(marker_table, file.path(P3, "10_marker_gene_detection.csv"), row.names = FALSE)
write.csv(panel_summary, file.path(P3, "11_marker_panel_coverage.csv"), row.names = FALSE)
write.csv(data.frame(gene = sort(shared_all)),
          file.path(P3, "12_genes_common_to_all_samples.csv"), row.names = FALSE)
write.csv(data.frame(gene = sort(union_all)),
          file.path(P3, "13_genes_in_any_sample.csv"), row.names = FALSE)

report <- c(
  "GSE150290 | Gate A Phase 3 | Gene feature compatibility audit",
  paste("Samples:", nrow(qc)),
  paste("Candidate independent patient IDs:", length(unique(qc$patient))),
  paste("Gene universe sizes:", paste(sort(unique(coverage$n_genes)), collapse = ", ")),
  paste("Intersection across all 52 samples:", length(shared_all)),
  paste("Union across all 52 samples:", length(union_all)),
  paste("Overlap between Pat01-A and Pat25-A:", length(intersect(ref_other, ref_early))),
  paste("All 47 Pat01-24 equal Pat01-A gene set:",
        all(coverage$identical_set_to_Pat01[grepl("^Pat(0[1-9]|1[0-9]|2[0-4])-", coverage$sample)])),
  paste("All 5 Pat25-29 equal Pat25-A gene set:",
        all(coverage$identical_set_to_Pat25[grepl("^Pat2[5-9]-", coverage$sample)])),
  "Do NOT interpret marker detection as proof of cell identity.",
  "Do NOT merge data until original sample pathology and feature processing are checked."
)
writeLines(report, file.path(P3, "14_phase3_summary.txt"))
cat("\n", paste(report, collapse = "\n"), "\n", sep = "")
cat("\nCandidate NAG/CAG marker panel coverage:\n")
print(subset(panel_summary,
             group_candidate %in% c("NAG_candidate", "CAG_candidate")),
      row.names = FALSE)
cat("\nOutputs: ", P3, "\n", sep = "")

}

# ---- legacy_02_Fig1_Donor_Exploration : Historical donor inference; DO NOT use as verified patient analysis ----
legacy_02_Fig1_Donor_Exploration <- function() {
# ============================================================================
# PLOS ONE | R1-4 / R2-1 | Fig1 patient/donor-level composition
# Runs on Windows, reads only current final 81,779-cell atlas metadata,
# and creates donor-level Supplementary XLSX + PDF/SVG + CSV audit outputs.
# No Seurat reclustering, no overwriting Fig1/Fig2/Fig3, no cell-level inference.
# ============================================================================

ROOT <- paste0(
  "E:/download/GSE150290_audit/GSE150290_fullcount_integration/",
  "Phase5_early_cohort/Phase7_repaired_and_integrated/",
  "Fig1_Original81779_Revised"
)
OUT <- "E:/download/GSE150290_audit/Fig1_Donor_Immune_20261021"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("openxlsx", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package unavailable: ", pkg)
}

csv <- function(x, name) utils::write.csv(
  x, file.path(OUT, name), row.names = FALSE, na = ""
)
error_out <- function(message, detail = NULL) {
  if (!is.null(detail)) csv(detail, "00_Donor_ID_REVIEW_BEFORE_ANALYSIS.csv")
  stop(message, call. = FALSE)
}

# 1. Prefer small existing cell annotation CSV to avoid loading 635-MB RDS
small <- file.path(ROOT, "All_Cell_Annotations.csv")
full <- file.path(ROOT, "Original81779_Fig1_Annotated_FINAL.rds")
needed <- c("sample_uid", "stage", "celltype_fig")
meta <- NULL
source_name <- NA_character_

if (file.exists(small)) {
  test <- tryCatch(utils::read.csv(small, stringsAsFactors = FALSE,
                                   check.names = FALSE), error = function(e) NULL)
  if (!is.null(test) && all(needed %in% names(test)) && nrow(test) == 81779L) {
    meta <- test
    source_name <- "All_Cell_Annotations.csv"
  }
  rm(test)
}
if (is.null(meta)) {
  if (!file.exists(full)) stop("Missing final atlas RDS: ", full)
  message("Loading 81,779-cell final Seurat object (metadata only needed) ...")
  obj <- readRDS(full)
  stopifnot(ncol(obj) == 81779L)
  meta <- obj@meta.data
  source_name <- "Original81779_Fig1_Annotated_FINAL.rds"
  rm(obj); gc(verbose = FALSE)
}
if (!all(needed %in% names(meta))) {
  error_out(paste("Missing metadata columns:",
                  paste(setdiff(needed, names(meta)), collapse = ", ")))
}
meta$sample_uid <- trimws(as.character(meta$sample_uid))
meta$stage <- toupper(trimws(as.character(meta$stage)))
meta$celltype <- as.character(meta$celltype_fig)
if (nrow(meta) != 81779L || anyNA(meta[, c("sample_uid", "stage", "celltype")])) {
  error_out("Cell count or required annotations invalid")
}

# 2. Determine donor IDs only from physical sample names, not provisional patient_id
# e.g. 5866_t1/t2 and 5866_n1/n2 => donor=5866; NAG1=>NAG1.
# If dataset-prefixed sample IDs exist, strip the :: prefix for donor parsing.
sample_map <- unique(meta[, c("sample_uid", "stage")])
if (anyDuplicated(sample_map$sample_uid)) {
  error_out("A sample_uid maps to multiple stages", sample_map)
}
sample_map$sample_name <- sub("^.*::", "", sample_map$sample_uid)
sample_map$donor_id <- ifelse(
  grepl("^[0-9]+_[nt][0-9]+$", sample_map$sample_name, ignore.case = TRUE),
  sub("_.*$", "", sample_map$sample_name),
  sample_map$sample_name
)
# Use dataset namespace when available (only if an ID collision is actually found)
if ("dataset" %in% names(meta)) {
  ds <- unique(meta[, c("sample_uid", "dataset")])
  if (anyDuplicated(ds$sample_uid)) error_out("Dataset ambiguous per sample", ds)
  sample_map$dataset <- as.character(ds$dataset[match(sample_map$sample_uid, ds$sample_uid)])
} else sample_map$dataset <- NA_character_
if ("patient_id" %in% names(meta)) {
  ps <- unique(meta[, c("sample_uid", "patient_id")])
  if (anyDuplicated(ps$sample_uid)) error_out("patient_id inconsistent within a sample", ps)
  sample_map$original_patient_id <- as.character(
    ps$patient_id[match(sample_map$sample_uid, ps$sample_uid)]
  )
} else sample_map$original_patient_id <- NA_character_

# Prevent a donor alias from accidentally merging samples across datasets
ds_valid <- subset(sample_map, !is.na(dataset) & dataset != "")
collision <- if (nrow(ds_valid)) {
  aggregate(dataset ~ donor_id, data = ds_valid,
            FUN = function(x) length(unique(x)))
} else data.frame(donor_id = character(), dataset = integer())
if (nrow(collision) && any(collision$dataset > 1L)) {
  error_out("Donor ID shared across datasets: manual audit required", sample_map)
}

# 3. Hard QC gates: verified 32 tissue samples and stage sample/donor counts
stages <- c("NAG", "CAG", "IM", "GC", "NAT")
observed <- data.frame(
  stage = stages,
  n_tissue_samples = vapply(stages, function(s) sum(sample_map$stage == s), integer(1)),
  n_donors = vapply(stages, function(s)
    length(unique(sample_map$donor_id[sample_map$stage == s])), integer(1))
)
expected_samples <- c(3L, 3L, 7L, 9L, 10L)
expected_donors <- c(3L, 3L, 7L, 7L, 8L)  # NAT=8 based on listed 10 adjacent tissues
csv(sample_map, "00_Donor_Sample_Mapping_AUDIT.csv")
csv(observed, "00_Stage_Unit_Counts_AUDIT.csv")
print(observed, row.names = FALSE)
if (nrow(sample_map) != 32L ||
    !all(as.integer(observed$n_tissue_samples) == expected_samples) ||
    !all(as.integer(observed$n_donors) == expected_donors) ||
    !setequal(unique(meta$stage), stages)) {
  error_out("STOP: samples/donors differ from known evidence; inspect 00_Donor... files", sample_map)
}

# 4. Explicitly preserve the 12-way celltype annotations; groups are not reclassified
alltypes <- sort(unique(meta$celltype))
required_types <- c("Macrophages", "T_cells", "B_cells", "NK_cells")
if (!all(required_types %in% alltypes) || length(alltypes) != 12L) {
  error_out(paste0("Expected Fig1 labels not found: ",
                   paste(setdiff(required_types, alltypes), collapse = ", ")),
            data.frame(actual_cell_types = alltypes))
}
meta$donor_id <- sample_map$donor_id[match(meta$sample_uid, sample_map$sample_uid)]
if (anyNA(meta$donor_id)) error_out("Some cells lack validated donor IDs")

# Full sample-level denominator includes all 12 types, not only immune types
sample_counts <- as.data.frame(table(
  sample_uid = factor(meta$sample_uid, levels = sample_map$sample_uid),
  celltype = factor(meta$celltype, levels = alltypes)
), stringsAsFactors = FALSE)
names(sample_counts)[3] <- "n_cells"
sample_counts$n_cells <- as.integer(sample_counts$n_cells)
sample_counts$stage <- sample_map$stage[match(sample_counts$sample_uid, sample_map$sample_uid)]
sample_counts$donor_id <- sample_map$donor_id[match(sample_counts$sample_uid, sample_map$sample_uid)]
sample_counts$sample_total <- ave(sample_counts$n_cells, sample_counts$sample_uid, FUN = sum)
sample_counts$proportion <- sample_counts$n_cells / sample_counts$sample_total

# Main analysis pools counts within stage/donor BEFORE dividing (avoids pseudoreplication)
donor_counts <- stats::aggregate(
  n_cells ~ stage + donor_id + celltype, sample_counts, FUN = sum
)
donor_total <- stats::aggregate(n_cells ~ stage + donor_id, donor_counts, FUN = sum)
names(donor_total)[3] <- "donor_total"
donor_counts <- merge(donor_counts, donor_total,
                      by = c("stage", "donor_id"), all.x = TRUE)
donor_counts$proportion <- donor_counts$n_cells / donor_counts$donor_total
sample_n <- stats::aggregate(sample_uid ~ stage + donor_id, sample_map,
                              FUN = function(x) length(unique(x)))
names(sample_n)[3] <- "n_tissue_samples"
donor_counts <- merge(donor_counts, sample_n,
                      by = c("stage", "donor_id"), all.x = TRUE)

# Sensitivity: equally weighted mean of constituent tissue sample proportions
sample_equal <- stats::aggregate(
  proportion ~ stage + donor_id + celltype, sample_counts, FUN = mean
)
names(sample_equal)[4] <- "mean_sample_proportion"
donor_counts <- merge(donor_counts, sample_equal,
                      by = c("stage", "donor_id", "celltype"), all.x = TRUE)
donor_counts$pool_minus_equal <- donor_counts$proportion - donor_counts$mean_sample_proportion

# 5. Four types used in comparisons; zero cells are retained as zero proportions
main <- donor_counts[donor_counts$celltype %in% required_types, ]
main$stage <- factor(main$stage, levels = stages)
main$celltype <- factor(main$celltype, levels = required_types)

f_summary <- function(v) c(
  n = length(v), mean = mean(v), sd = if(length(v)>1) stats::sd(v) else NA_real_,
  median = stats::median(v), Q1 = unname(stats::quantile(v, 0.25)),
  Q3 = unname(stats::quantile(v, 0.75)), min = min(v), max = max(v)
)
stage_summary <- do.call(rbind, lapply(required_types, function(typ) {
  do.call(rbind, lapply(stages, function(st) {
    d <- subset(main, celltype == typ & stage == st)
    cbind(data.frame(stage = st, celltype = typ),
          as.data.frame(as.list(f_summary(d$proportion))))
  }))
}))
rownames(stage_summary) <- NULL
stage_summary$unit <- "donor within stage; pooled cell-count fraction"

# 6. Exploratory unpaired stage comparisons, NOT single-cell hypothesis testing.
# NAT excluded from ordered comparisons. If donors overlap across stages,
# never silently treat them as independent.
ordered <- c("NAG", "CAG", "IM", "GC")
pairs <- combn(ordered, 2, simplify = FALSE)
comparisons <- do.call(rbind, lapply(required_types, function(typ) {
  do.call(rbind, lapply(pairs, function(pair) {
    a <- subset(main, celltype == typ & stage == pair[1])
    b <- subset(main, celltype == typ & stage == pair[2])
    overlap <- intersect(as.character(a$donor_id), as.character(b$donor_id))
    p <- if(length(overlap)) NA_real_ else
      suppressWarnings(stats::wilcox.test(a$proportion, b$proportion,
                                           exact = FALSE)$p.value)
    data.frame(celltype = typ, stage_A = pair[1], stage_B = pair[2],
      n_A = nrow(a), n_B = nrow(b),
      mean_A = mean(a$proportion), mean_B = mean(b$proportion),
      median_A = median(a$proportion), median_B = median(b$proportion),
      mean_difference_B_minus_A = mean(b$proportion) - mean(a$proportion),
      median_difference_B_minus_A = median(b$proportion) - median(a$proportion),
      n_shared_donors = length(overlap), p_wilcoxon_exploratory = p,
      note = if(length(overlap)) "NOT TESTED: shared donor across stages" else
        "Exploratory; limited n; independent donor units"
    )
  }))
}))
comparisons$q_BH_all_fourtypes_24tests <- stats::p.adjust(
  comparisons$p_wilcoxon_exploratory, method = "BH"
)

# 7. GC and NAT have donor overlap; paired difference shown descriptively
wGC <- subset(main, stage == "GC", select = c("donor_id", "celltype", "proportion"))
wNAT <- subset(main, stage == "NAT", select = c("donor_id", "celltype", "proportion"))
paired <- merge(wGC, wNAT, by = c("donor_id", "celltype"),
                suffixes = c("_GC", "_NAT"))
paired$delta_GC_minus_NAT <- paired$proportion_GC - paired$proportion_NAT

# 8. Save complete, styled XLSX
wb <- openxlsx::createWorkbook()
styles <- list(
  hdr = openxlsx::createStyle(fgFill = "#17324D", fontColour = "#FFFFFF",
                               textDecoration = "bold", halign = "center",
                               valign = "center", wrapText = TRUE),
  num = openxlsx::createStyle(numFmt = "0.0000"),
  pct = openxlsx::createStyle(numFmt = "0.00%")
)
addsheet <- function(name, df, percent_cols = character()) {
  openxlsx::addWorksheet(wb, name, gridLines = FALSE)
  openxlsx::writeData(wb, name, df, headerStyle = styles$hdr, withFilter = TRUE)
  openxlsx::freezePane(wb, name, firstRow = TRUE)
  # Capped widths keep the supplementary workbook printable and readable
  wide <- ifelse(nchar(names(df)) > 23, 26, 17)
  wide[grepl("note|Details|sample_uid|patient|donor_id|unit", names(df), ignore.case=TRUE)] <- 31
  openxlsx::setColWidths(wb, name, cols = seq_along(df), widths = wide)
  # Keep text columns from becoming excessively wide
  openxlsx::setRowHeights(wb, name, rows = 1, heights = 42)
  if(nrow(df) && length(percent_cols)) {
    for(nm in intersect(percent_cols, names(df))) {
      j <- match(nm, names(df))
      openxlsx::addStyle(wb, name, styles$pct, rows = 2:(nrow(df)+1),
                          cols = j, gridExpand = TRUE, stack = TRUE)
    }
  }
}
notes <- data.frame(
  Field = c("Atlas", "Input", "Main unit", "Donor ID", "GC replicate handling",
            "NAT", "Denominator", "Primary proportion", "Sensitivity measure",
            "Hypothesis testing", "Multiplicity", "Interpretation", "Important limitation"),
  Details = c("81,779 cells / 32 tissue samples", source_name,
              "Independent donor within disease stage; not individual cells",
              "Numeric sample prefix for t/n biopsies; singleton sample ID otherwise",
              "Pool counts from 5866_t1/t2 and 5931_t1/t2 to donor level",
              "Adjacent-tissue reference, not in ordered disease-stage tests",
              "All annotated cells (12 types), NOT selected immune-only sum",
              "Sum cell counts of same donor/stage/type divided by all cells from donor/stage",
              "Mean of per-tissue proportions for donors with multiple samples",
              "Exploratory two-sided Wilcoxon; shared donor pairs not tested as independent",
              "BH correction across up to 24 four-lineage NAG/CAG/IM/GC comparisons",
              "Descriptive/association, not proof of abundance changes or causal development",
              "NAG and CAG n=3; compositional fractions; dataset confounding; no absolute counts")
)
addsheet("README_Methods", notes)
addsheet("Donor_Sample_Map", sample_map)
addsheet("Stage_Units", observed)
addsheet("Sample_Counts_All12", sample_counts, "proportion")
addsheet("Donor_Counts_All12", donor_counts,
         c("proportion", "mean_sample_proportion", "pool_minus_equal"))
addsheet("Donor_Immune_4", main,
         c("proportion", "mean_sample_proportion", "pool_minus_equal"))
addsheet("Stage_Summary_4", stage_summary,
         c("mean", "sd", "median", "Q1", "Q3", "min", "max"))
addsheet("Pairwise_Exploratory", comparisons,
         c("mean_A", "mean_B", "median_A", "median_B",
           "mean_difference_B_minus_A", "median_difference_B_minus_A"))
addsheet("GC_NAT_Paired", paired,
         c("proportion_GC", "proportion_NAT", "delta_GC_minus_NAT"))

xlsx <- file.path(OUT, "Supplementary_DonorLevel_Immune_Composition.xlsx")
openxlsx::saveWorkbook(wb, xlsx, overwrite = TRUE)

# 9. Independent donor dot+box plots; NAT is a visually separate reference
set.seed(20261021)
plotdat <- main
plotdat$stage <- factor(plotdat$stage, levels = stages)
plotdat$celltype <- factor(plotdat$celltype,
                           levels = c("Macrophages", "T_cells", "B_cells", "NK_cells"))
cols <- c(NAG = "#5C7B91", CAG = "#D8A35D", IM = "#65A694",
          GC = "#C66B76", NAT = "#8A8F98")
p <- ggplot2::ggplot(plotdat, ggplot2::aes(stage, proportion, colour = stage)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = .48, colour = "#505050",
                        fill = "white", linewidth = .4) +
  ggplot2::geom_jitter(width = .11, height = 0, size = 1.9, alpha = .95) +
  ggplot2::facet_wrap(~ celltype, ncol = 2) +
  ggplot2::scale_colour_manual(values = cols, guide = "none") +
  ggplot2::scale_y_continuous(limits = c(0, 1),
    labels = function(x) paste0(round(100*x), "%")) +
  ggplot2::labs(x = "Disease stage / reference tissue",
                y = "Proportion among all captured cells",
                title = "Donor-level immune-cell composition",
                subtitle = "One point per donor per stage; GC repeated biopsies pooled") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                 plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(OUT, "Supplementary_Immune_DonorLevel.pdf"),
                p, width = 10, height = 7)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggplot2::ggsave(file.path(OUT, "Supplementary_Immune_DonorLevel.svg"),
                  p, width = 10, height = 7)
}

# 10. Audit all n and cell totals
stopifnot(sum(sample_counts$n_cells) == 81779L,
          sum(donor_counts$n_cells) == 81779L,
          nrow(main) == sum(observed$n_donors) * 4L)
csv(donor_counts, "Donor_Counts_All12.csv")
csv(main, "Donor_Immune_4.csv")
csv(comparisons, "Pairwise_Exploratory.csv")
cat("\nPASS: 81,779 cells; 32 tissue samples; donor mapping checked.\n")
cat("XLSX:", xlsx, "\nFIGURE:", file.path(OUT, "Supplementary_Immune_DonorLevel.pdf"), "\n")
cat("Note: Exploratory comparisons are not proof of biological replication across datasets.\n")

}

# ---- legacy_03_Fig2_CellChat_HistoricalPlot : Historical Fig2 plotting; OLD D: paths ----
legacy_03_Fig2_CellChat_HistoricalPlot <- function() {
# Fig2 v2 — Compact SCI-style composite (PLOS ONE revision)
# Real data from 4 existing CellChat .RData files. No synthetic values.
# Improvements: aligned four-matrix row; pooled nonzero signal bins;
# equal-height B and C panels; compact labeling and minimal wasted space.
# Adjust root/out_dir in CONFIG. Requires CellChat, ggplot2, patchwork.

# ==============================================================
# PLOS ONE revision (PONE-D-26-36895)
# Figure 2: Comparative intercellular communication across stages
# Reads existing CellChat .RData files; does NOT alter their objects.
# Output A: 4-stage directional communication heatmaps
# Output B: cell-type outgoing/incoming share of network strength
# Output C: directed macrophage-associated ligand-receptor dot plot
# IMPORTANT: descriptive pooled-stage analysis, not patient-level inference.
# ==============================================================

# -------- 0. Configure paths (edit only these lines) --------------
root <- "D:/陈建/本科毕设/10.22号 结果/cellchat"
out_dir <- "D:/陈建/本科毕设/4.5 结果/PLOS_ONE_revision/Fig2_revision"
files <- c(
  NAG = file.path(root, "NAG", "NAG.RData"),
  CAG = file.path(root, "CAG", "CAG.RData"),
  IM  = file.path(root, "IM", "IM.RData"),
  GC  = file.path(root, "Tumor", "Tumor.RData")
)
# If the filenames differ, edit the entries in 'files'.

missing <- files[!file.exists(files)]
if (length(missing)) stop("Missing .RData files:\n", paste(missing, collapse = "\n"))
for (pkg in c("CellChat", "ggplot2", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not installed: ", pkg,
         ". Install dependencies first; see Fig2_README.txt")
  }
}
suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stages <- c("NAG", "CAG", "IM", "GC")
stage_col <- c(NAG = "#477B5D", CAG = "#91A567",
               IM = "#C48C51", GC = "#B8564D")

# Different spellings from saved CellChat objects are standardized for DISPLAY
# only. Stop on ambiguous collisions rather than silently merging populations.
canonical_cell <- function(x) {
  raw <- trimws(gsub("_", " ", as.character(x), fixed = TRUE))
  key <- tolower(gsub("[^[:alnum:]]", "", raw))
  dict <- c(
    macrophages="Macrophages", macrophage="Macrophages", mac="Macrophages",
    epithelialcells="Epithelial cells", epithelial="Epithelial cells", epi="Epithelial cells",
    fibroblasts="Fibroblasts", fibroblast="Fibroblasts",
    tcells="T cells", tcell="T cells", t="T cells",
    nkcells="NK cells", nkcell="NK cells", nk="NK cells",
    bcells="B cells", bcell="B cells", b="B cells",
    plasmacells="Plasma cells", plasmacell="Plasma cells", plasma="Plasma cells",
    mastcells="Mast cells", mastcell="Mast cells", mast="Mast cells",
    endothelialcells="Endothelial cells", endothelial="Endothelial cells",
    vascularsmoothmusclecells="Vascular smooth muscle cells"
  )
  out <- raw
  hit <- key %in% names(dict)
  out[hit] <- unname(dict[key[hit]])
  out
}

# Load into an isolated environment; do not overwrite NAG, Tumor, a1, etc.
read_cellchat <- function(path, stage) {
  env <- new.env(parent = .GlobalEnv)
  restored <- load(path, envir = env)
  is_cc <- vapply(restored, function(nm) {
    inherits(get(nm, envir = env), "CellChat")
  }, logical(1))
  opts <- restored[is_cc]
  if (!length(opts)) {
    stop("No CellChat object in ", path,
         ". Saved objects: ", paste(restored, collapse = ", "))
  }
  if (length(opts) > 1L) {
    if ("cellchat" %in% opts) opts <- "cellchat" else {
      stop("Multiple CellChat objects in ", path, ": ",
           paste(opts, collapse = ", "), ". Please select explicitly.")
    }
  }
  message(stage, ": recovered object '", opts[1], "' from ", basename(path))
  get(opts[1], envir = env)
}

get_stage_result <- function(path, stage) {
  cc <- read_cellchat(path, stage)
  m <- cc@net$weight
  if (is.null(m) || !is.matrix(m)) {
    if (is.null(cc@net$prob) || is.null(cc@net$pval)) {
      stop(stage, ": missing both network weight matrix and probability arrays")
    }
    cc <- CellChat::aggregateNet(cc, thresh = 0.05)
    m <- cc@net$weight
  }
  if (nrow(m) != ncol(m)) stop(stage, ": non-square network matrix")
  if (is.null(rownames(m)) || is.null(colnames(m))) {
    stop(stage, ": the weight matrix has no cell-type names")
  }
  rownames(m) <- canonical_cell(rownames(m))
  colnames(m) <- canonical_cell(colnames(m))
  if (anyDuplicated(rownames(m)) || anyDuplicated(colnames(m))) {
    stop(stage, ": duplicate cell-type labels after spelling standardization")
  }
  if (!setequal(rownames(m), colnames(m))) {
    stop(stage, ": source and target cell types differ")
  }
  if (any(!is.finite(m)) || any(m < 0)) {
    stop(stage, ": invalid network strengths")
  }
  m <- m[rownames(m), rownames(m), drop = FALSE]

  # CellChat bootstrap p-values are within-stage, NOT between-stage tests.
  ed <- tryCatch(
    CellChat::subsetCommunication(cc, slot.name = "net", thresh = 0.05),
    error = function(e) {
      stop(stage, ": could not extract ligand-receptor results: ",
           conditionMessage(e))
    }
  )
  if (is.null(ed)) ed <- data.frame()
  if (nrow(ed) > 0L) {
    needed <- c("source", "target", "prob")
    if (!all(needed %in% names(ed)))
      stop(stage, ": LR output is missing ",
           paste(setdiff(needed, names(ed)), collapse = ", "))
    ed$source <- canonical_cell(ed$source)
    ed$target <- canonical_cell(ed$target)
    ed$stage <- stage
    ed$prob <- as.numeric(ed$prob)
    if ("interaction_name_2" %in% names(ed)) {
      ed$LR_label <- as.character(ed$interaction_name_2)
    } else if (all(c("ligand", "receptor") %in% names(ed))) {
      ed$LR_label <- paste(ed$ligand, ed$receptor, sep = " - ")
    } else if ("interaction_name" %in% names(ed)) {
      ed$LR_label <- as.character(ed$interaction_name)
    } else stop(stage, ": no ligand-receptor label columns found")
    ed$pathway_name <- if ("pathway_name" %in% names(ed))
      as.character(ed$pathway_name) else ""
  }

  # Record the reference LR database where available.
  dbnames <- tryCatch(as.character(cc@DB$interaction$interaction_name),
                      error = function(e) character())
  param_text <- paste(capture.output(str(cc@options, max.level = 2)),
                      collapse = "\n")
  data.frame_summary <- data.frame(
    stage = stage,
    n_cell_types = nrow(m),
    n_cells_in_saved_object = length(cc@idents),
    total_network_strength = sum(m),
    n_significant_LR_rows = nrow(ed),
    stringsAsFactors = FALSE
  )
  list(weight = m, edges = ed, summary = data.frame_summary,
       dbnames = sort(unique(dbnames)), params = param_text)
}

message("Reading the four saved stage-level CellChat objects...")
result <- lapply(stages, function(s) get_stage_result(files[[s]], s))
names(result) <- stages
summary_df <- do.call(rbind, lapply(result, `[[`, "summary"))
write.csv(summary_df, file.path(out_dir, "Fig2_stage_network_audit.csv"),
          row.names = FALSE)
writeLines(unlist(lapply(stages, function(s) {
  c(paste0("==== ", s, " ===="), result[[s]]$params, "")
})), file.path(out_dir, "Fig2_saved_CellChat_options.txt"))

# Checks before cross-stage visualization.
ref_db <- result[[1]]$dbnames
for (s in stages[-1]) {
  if (length(ref_db) && length(result[[s]]$dbnames) &&
      !identical(ref_db, result[[s]]$dbnames)) {
    warning(s, ": CellChat LR database entries differ from NAG; ",
            "interpret cross-stage comparisons cautiously")
  }
}

cell_priority <- c("Epithelial cells", "T cells", "NK cells", "B cells",
                   "Plasma cells", "Macrophages", "Mast cells",
                   "Fibroblasts", "Endothelial cells",
                   "Vascular smooth muscle cells")
all_cells <- unique(unlist(lapply(result, function(z) rownames(z$weight))))
all_cells <- c(intersect(cell_priority, all_cells),
               sort(setdiff(all_cells, cell_priority)))


# ================================================================
# v2 | Balanced journal figure. Uses saved inference; no CellChat
# re-inference. Data values are not altered.
# ================================================================

# Set global publication font to an available serif font.
# This can be changed to "Times New Roman" if available on Windows.
FONT <- "serif"
FIG_WIDTH <- 15.6
FIG_HEIGHT <- 10.0
DPI <- 450

# ---------- Common labels and complete aligned stage matrices ----
short_names <- c(
  "Epithelial cells" = "Epi", "T cells" = "T", "NK cells" = "NK",
  "B cells" = "B", "Plasma cells" = "Plasma", "Macrophages" = "Mac",
  "Mast cells" = "Mast", "Fibroblasts" = "Fib",
  "Endothelial cells" = "Endo", "Vascular smooth muscle cells" = "VSMC"
)
short_cell <- function(x) {
  x <- as.character(x)
  out <- x
  hit <- x %in% names(short_names)
  out[hit] <- unname(short_names[x[hit]])
  out
}
short_levels <- short_cell(all_cells)
if (anyDuplicated(short_levels)) stop("Cell names collide after abbreviation")

mat_to_df <- function(z, stage) {
  m <- z$weight
  total <- sum(m)
  if (total <= 0) stop(stage, " has zero total strength")
  d <- as.data.frame(as.table(100 * m / total), stringsAsFactors=FALSE)
  colnames(d) <- c("source", "target", "strength_percent")
  d$stage <- stage
  d
}
net_df <- do.call(rbind, Map(mat_to_df, result, stages))
full <- expand.grid(stage=stages, source=all_cells, target=all_cells,
                    stringsAsFactors=FALSE)
net_df <- merge(full, net_df, by=c("stage", "source", "target"), all.x=TRUE,
                sort=FALSE)
net_df$stage <- factor(net_df$stage, levels=stages)
net_df$source_short <- factor(short_cell(net_df$source), levels=short_levels)
net_df$target_short <- factor(short_cell(net_df$target), levels=rev(short_levels))
if (any(net_df$strength_percent < 0, na.rm=TRUE)) stop("Negative network weight")

# IMPORTANT: one pooled, stage-shared scale. 0 stays 0 (no connection).
# Positive weights are ranked into shared quintiles, rather than
# being washed out by one Epi->Epi high-strength connection.
positive <- net_df$strength_percent[is.finite(net_df$strength_percent) &
                                      net_df$strength_percent > 0]
if (length(positive) < 5) stop("Too few nonzero interactions for a figure")
q <- as.numeric(stats::quantile(positive, c(.2, .4, .6, .8), na.rm=TRUE))
net_df$strength_bin <- vapply(net_df$strength_percent, function(x) {
  if (is.na(x)) return("Unavailable")
  if (x == 0) return("0")
  paste0("Q", min(5L, 1L+findInterval(x, q)))
}, character(1))
bin_levels <- c("0", paste0("Q",1:5), "Unavailable")
net_df$strength_bin <- factor(net_df$strength_bin, levels=bin_levels)

palette_A <- c(
  "0" = "#F3F5F6", "Q1" = "#D7E5EB", "Q2" = "#ACCBD7",
  "Q3" = "#79AABD", "Q4" = "#397D9C", "Q5" = "#164867",
  "Unavailable" = "#D8D8D8"
)
write.csv(data.frame(
  interval=c("0", "Q1: lowest 20% nonzero", "Q2: 20%-40%",
             "Q3: 40%-60%", "Q4: 60%-80%", "Q5: top 20%", "Unavailable"),
  pooled_strength_percent_reference=c(0, q[1], q[2], q[3], q[4], max(positive), NA),
  stringsAsFactors=FALSE
), file.path(out_dir, "Fig2A_pooled_strength_scale_v2.csv"), row.names=FALSE)
write.csv(net_df, file.path(out_dir, "Fig2A_stage_network_values_v2.csv"), row.names=FALSE)

# Use a single faceted ggplot; guarantees all four matrices have
# identical dimensions, axis scales and consistent row/column order.
base_theme <- theme_classic(base_family=FONT, base_size=10) +
  theme(
    axis.line=element_blank(),
    axis.ticks=element_blank(),
    panel.grid=element_blank(),
    strip.background=element_rect(fill="#F3F5F6", colour="#D9DDE0", linewidth=.30),
    strip.text=element_text(face="bold", size=11, margin=margin(4,2,4,2)),
    legend.title=element_text(size=9.5, face="bold"),
    legend.text=element_text(size=8.5),
    plot.title=element_text(face="bold", size=13, margin=margin(0,0,6,0)),
    plot.subtitle=element_blank(),
    plot.margin=margin(3,4,2,3,unit="pt")
  )

pA <- ggplot(net_df, aes(x=source_short, y=target_short, fill=strength_bin)) +
  geom_tile(colour="white", linewidth=.45, width=.98, height=.98) +
  facet_grid(. ~ stage, drop=FALSE) +
  scale_fill_manual(values=palette_A, drop=FALSE,
                    breaks=c("0","Q1","Q2","Q3","Q4","Q5"),
                    labels=c("0", "Q1", "Q2", "Q3", "Q4", "Q5"),
                    name="Shared edge\nstrength bin") +
  scale_x_discrete(drop=FALSE, expand=expansion(mult=0)) +
  scale_y_discrete(drop=FALSE, expand=expansion(mult=0)) +
  coord_fixed(ratio=1, clip="off") +
  guides(fill=guide_legend(override.aes=list(colour=NA), ncol=1,
                           keyheight=grid::unit(0.32,"cm"))) +
  labs(title="A   Cell-to-cell communication patterns",
       x="Sender", y="Receiver") +
  base_theme +
  theme(axis.text.x=element_text(size=8.5, angle=45, hjust=1, vjust=1),
        axis.text.y=element_text(size=8.5),
        axis.title=element_text(size=9.5),
        panel.spacing.x=grid::unit(0.18, "cm"),
        legend.position="right", legend.box.margin=margin(0,0,0,0))

# --------- Panel B: comparative sender/receiver contribution ------
centrality_df <- do.call(rbind, lapply(stages, function(s) {
  m <- result[[s]]$weight
  total <- sum(m)
  if (total <= 0) stop("Zero network strength in ", s)
  data.frame(stage=s,
             celltype=rep(rownames(m),2),
             direction=rep(c("Sender", "Receiver"), each=nrow(m)),
             strength_share_percent=100*c(rowSums(m),colSums(m))/total,
             stringsAsFactors=FALSE)
}))
centrality_df$stage <- factor(centrality_df$stage,levels=stages)
centrality_df$cell_short <- factor(short_cell(centrality_df$celltype), levels=rev(short_levels))
centrality_df$direction <- factor(centrality_df$direction,levels=c("Sender","Receiver"))
centrality_df$highlight <- ifelse(centrality_df$celltype=="Macrophages", "Mac", "Other")
write.csv(centrality_df,file.path(out_dir,"Fig2B_sender_receiver_contributions_v2.csv"),row.names=FALSE)

pB <- ggplot(centrality_df, aes(x=stage,y=cell_short)) +
  geom_tile(fill="#FAFBFC",color="#E7EAEC",linewidth=.38,width=.97,height=.97) +
  geom_point(aes(size=strength_share_percent,fill=highlight),
             shape=21, colour="#66737B",stroke=.28) +
  facet_grid(.~direction, drop=FALSE) +
  scale_size_area(max_size=7.5, breaks=c(5,15,30,50),
                  name="Network share (%)") +
  scale_fill_manual(values=c(Mac="#C26755",Other="#A5B4BE"),guide="none") +
  scale_x_discrete(drop=FALSE,expand=expansion(mult=.02)) +
  scale_y_discrete(drop=FALSE,expand=expansion(mult=.02)) +
  labs(title="B   Cell-type communication contributions", x=NULL,y=NULL) +
  guides(size=guide_legend(override.aes=list(fill="#8EA2AC"),nrow=1)) +
  base_theme +
  theme(axis.text=element_text(size=8.7),
        axis.text.y=element_text(size=8.7),
        strip.text=element_text(size=10),
        panel.spacing.x=grid::unit(0.12,"cm"),
        legend.position="bottom",legend.direction="horizontal",
        legend.box.margin=margin(0,0,0,0))

# --------- Panel C: balanced, directed macrophage signaling -------
ed_list <- lapply(stages,function(s) result[[s]]$edges)
ed_list <- Filter(function(z) nrow(z)>0,ed_list)
if (!length(ed_list)) stop("Saved objects have no extracted LR edges")
ed <- do.call(rbind,ed_list)
ed <- ed[ed$source=="Macrophages" | ed$target=="Macrophages",,drop=FALSE]
if (!nrow(ed)) stop("No macrophage-associated LR detected")
ed$edge_label <- paste0(short_cell(ed$source), " -> ", short_cell(ed$target),
                         "  |  ", ed$LR_label)
# preselected signatures are represented if detected, but no pathway is forced.
ed$family <- ifelse(grepl("MIF",ed$pathway_name,ignore.case=TRUE) |
                       grepl("MIF",ed$LR_label,ignore.case=TRUE), "MIF",
                    ifelse(grepl("TNFSF13B|BAFF",ed$pathway_name,ignore.case=TRUE) |
                           grepl("TNFSF13B|BAFF",ed$LR_label,ignore.case=TRUE),
                           "BAFF", "Other"))
# prioritize reproducibility in pooled stages, not simply top one dominant pathway.
by_edge <- aggregate(prob ~ edge_label + family,ed,sum)
colnames(by_edge)[3] <- "total_probability"
ns <- aggregate(stage ~ edge_label, unique(ed[,c("edge_label","stage")]),
                function(z)length(unique(z)))
colnames(ns)[2] <- "n_stages"
by_edge <- merge(by_edge,ns,by="edge_label",all.x=TRUE)
by_edge <- by_edge[order(-by_edge$n_stages,-by_edge$total_probability),,drop=FALSE]
# up to 4 MIF, 2 BAFF, 4 other signatures. Unused slots are filled by remaining signals.
chosen <- character()
for (family in c("MIF","BAFF","Other")) {
  cap <- c(MIF=4L,BAFF=2L,Other=4L)[[family]]
  candidates <- by_edge$edge_label[by_edge$family==family]
  chosen <- c(chosen,head(setdiff(candidates,chosen),cap))
}
chosen <- unique(c(chosen,head(setdiff(by_edge$edge_label,chosen),10L-length(chosen))))
chosen <- head(chosen,10L)
if (!length(chosen)) stop("No candidate ligand receptor interactions")

ed_sel <- ed[ed$edge_label %in% chosen,,drop=FALSE]
# Duplicate stages/edges are aggregated for display, preserving all probabilities.
plot_ed <- aggregate(prob~stage+edge_label,data=ed_sel,sum)
full_ed <- expand.grid(stage=stages,edge_label=chosen,stringsAsFactors=FALSE)
plot_ed <- merge(full_ed,plot_ed,by=c("stage","edge_label"),all.x=TRUE,sort=FALSE)
plot_ed$stage <- factor(plot_ed$stage,levels=stages)
plot_ed$edge_label <- factor(plot_ed$edge_label,levels=rev(chosen))
write.csv(plot_ed,file.path(out_dir,"Fig2C_selected_macrophage_LR_v2.csv"),row.names=FALSE)

# Keep stage colors consistent with the user's other figures.
pC <- ggplot(plot_ed,aes(x=stage,y=edge_label)) +
  geom_tile(fill="#FAFBFC",colour="#E7EAEC",linewidth=.38,width=.97,height=.97) +
  geom_point(data=subset(plot_ed,!is.na(prob)),
             aes(size=prob,fill=stage),shape=21,color="#4A555B",stroke=.26) +
  scale_size_continuous(range=c(2.4,7),name="Inferred LR\nprobability") +
  scale_fill_manual(values=stage_col,guide="none") +
  scale_x_discrete(drop=FALSE,expand=expansion(mult=.025)) +
  scale_y_discrete(drop=FALSE,expand=expansion(mult=.025)) +
  labs(title="C   Macrophage-associated LR signaling",x=NULL,y=NULL) +
  guides(size=guide_legend(nrow=1,override.aes=list(fill="#618E91"))) +
  base_theme +
  theme(axis.text.x=element_text(size=9),axis.text.y=element_text(size=8),
        legend.position="bottom",legend.direction="horizontal",
        legend.box.margin=margin(0,0,0,0))

# ---------- Single compact, rigorously aligned 2-row layout ------
# A across full width; B and C share the same second-row height.
# No automatic letter tags which alter facet-panel measurements.
bottom <- (pB | pC) + patchwork::plot_layout(widths=c(1.05,1.37))
combined <- (pA / bottom) +
  patchwork::plot_layout(heights=c(.98,1.09)) +
  patchwork::plot_annotation(
    theme=theme(plot.background=element_rect(fill="white",colour=NA),
                plot.margin=margin(4,4,4,4,unit="pt")))

outpdf <- file.path(out_dir,"Fig2_revised_v2.pdf")
outtif <- file.path(out_dir,"Fig2_revised_v2.tif")
ggsave(outpdf,combined,width=FIG_WIDTH,height=FIG_HEIGHT,units="in",
       limitsize=FALSE,bg="white",device="pdf")
ggsave(outtif,combined,width=FIG_WIDTH,height=FIG_HEIGHT,units="in",
       dpi=DPI,compression="lzw",limitsize=FALSE,bg="white",device="tiff")
# Separate panels allow layout adjustment without distorting figure elements.
ggsave(file.path(out_dir,"Fig2A_v2.pdf"),pA,width=15.6,height=4.9,units="in",bg="white")
ggsave(file.path(out_dir,"Fig2B_v2.pdf"),pB,width=6.7,height=5.25,units="in",bg="white")
ggsave(file.path(out_dir,"Fig2C_v2.pdf"),pC,width=8.5,height=5.25,units="in",bg="white")

cat("\nSaved publication layout:\n",outpdf,"\n",outtif,"\n")
cat("A heatmap positive edges:",length(positive),
    "; shared nonzero quintiles (%):",paste(signif(q,3),collapse=", "),"\n")
cat("Selected unique LR pairs:",length(chosen),"\n")
cat("No patient-level p-values are inferred from these pooled stage-level networks.\n")

}

# ---- legacy_04_CellChat_CIB_SavedAudit : Saved-object audit ----
legacy_04_CellChat_CIB_SavedAudit <- function() {
# PLOS ONE R1-5 / R2-4: audit saved CellChat and existing CIBERSORT outputs.
# No CellChat re-fitting, no CIBERSORT re-running, no changes to Fig2/Fig3.
ROOT <- paste0("E:/download/GSE150290_audit/",
  "GSE150290_fullcount_integration/Phase5_early_cohort/",
  "Phase7_repaired_and_integrated/Fig1_Original81779_Revised")
OUT <- "E:/download/GSE150290_audit/CellChat_CIB_Audit_20261021"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("CellChat", "openxlsx")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("R package missing: ", pkg)
}
writecsv <- function(df, name) write.csv(df, file.path(OUT, name), row.names=FALSE, na="")
stages <- c("NAG", "CAG", "IM", "GC")
ccdir <- file.path(ROOT,"Fig2_CellChat_Revised")
ccpaths <- file.path(ccdir, paste0("CellChat_", stages, "_primary.rds"))
if (!all(file.exists(ccpaths))) {
  stop("Original CellChat stage objects missing; inspect: ", ccdir)
}
getp <- function(o,k) {
  v <- o@options$parameter[[k]]
  if(is.null(v) || !length(v)) NA_character_ else paste(v,collapse=";")
}
opts <- list(); centers <- list(); nets <- list()
for(i in seq_along(stages)) {
  st <- stages[i]
  message("Reading saved CellChat: ", st)
  obj <- readRDS(ccpaths[i])
  weight <- obj@net$weight
  if(is.null(weight) || length(dim(weight))!=2L || nrow(weight)!=ncol(weight)) {
    stop("Missing square @net$weight for ",st)
  }
  weight <- as.matrix(weight)
  if (!is.null(rownames(weight)) && !is.null(colnames(weight)) &&
      !identical(rownames(weight), colnames(weight))) {
    if (setequal(rownames(weight), colnames(weight))) {
      weight <- weight[, rownames(weight), drop = FALSE]
    } else stop("Sender/receiver labels incompatible in ", st)
  }
  weight[!is.finite(weight)] <- 0
  nc <- obj@net$count
  has_count <- !is.null(nc) && identical(dim(nc),dim(weight))
  if (!has_count) nc <- matrix(NA_real_,nrow(weight),ncol(weight))
  nc <- as.matrix(nc)
  ct <- rownames(weight)
  if(is.null(ct)) ct <- colnames(weight)
  if(is.null(ct)) {
    ct <- levels(obj@idents)
    if(length(ct)!=nrow(weight)) stop("Cannot identify network labels for ",st)
  }
  db_interactions <- tryCatch(nrow(obj@DB$interaction),error=function(e) NA_integer_)
  opts[[i]] <- data.frame(
    stage=st,n_cells=length(obj@idents),n_types=length(ct),
    type.mean=getp(obj,"type.mean"),trim=getp(obj,"trim"),
    raw.use=getp(obj,"raw.use"),population.size=getp(obj,"population.size"),
    nboot=getp(obj,"nboot"),Kh=getp(obj,"Kh"),
    seed.use=getp(obj,"seed.use"),
    min.cells="UNVERIFIED: NOT stored in @options$parameter",
    DB_release="UNVERIFIED: source code required",
    n_DB_interactions=db_interactions
  )
  outgoing <- rowSums(weight)
  incoming <- colSums(weight)
  nn <- nrow(weight)
  total_weight <- sum(weight)
  total_counts <- if(has_count) sum(nc, na.rm=TRUE) else NA_real_
  x <- data.frame(
    stage=st,celltype=ct,
    outgoing_weight=as.numeric(outgoing), incoming_weight=as.numeric(incoming),
    total_weight=as.numeric(outgoing+incoming),
    outgoing_count=if(has_count) as.numeric(rowSums(nc,na.rm=TRUE)) else rep(NA_real_,nn),
    incoming_count=if(has_count) as.numeric(colSums(nc,na.rm=TRUE)) else rep(NA_real_,nn),
    outgoing_rank=rank(-outgoing,ties.method="min"),
    incoming_rank=rank(-incoming,ties.method="min"),
    outgoing_fraction_of_stage=if(total_weight>0) outgoing/total_weight else NA_real_,
    incoming_fraction_of_stage=if(total_weight>0) incoming/total_weight else NA_real_,
    n_types_in_stage=nn,
    stringsAsFactors=FALSE
  )
  centers[[i]] <- x
  nets[[i]] <- data.frame(stage=st,n_cells=length(obj@idents),
     n_celltypes=nn,n_positive_edges=sum(weight>0),
     inferred_count_sum=total_counts,inferred_weight_sum=total_weight)
  rm(obj);gc(verbose=FALSE)
}
parameters <- do.call(rbind,opts)
centrality <- do.call(rbind,centers)
network <- do.call(rbind,nets)
writecsv(parameters,"01_CellChat_Saved_Parameters.csv")
writecsv(centrality,"02_CellChat_WithinStage_Strength_Rank.csv")
writecsv(network,"03_CellChat_Saved_Networks.csv")

# Previously generated descriptive comparison. These sources are NOT a fresh reanalysis.
d3 <- file.path(ROOT,"Fig3_ABC_Revised")
trends_file <- file.path(d3,"Fig3E_CrossPlatform_Trends.csv")
sample_file <- file.path(d3,"Fig3E_CrossPlatform_Sample_Data.csv")
trends <- if(file.exists(trends_file)) read.csv(trends_file,check.names=FALSE) else NULL
samples <- if(file.exists(sample_file)) read.csv(sample_file,check.names=FALSE) else NULL
if(is.null(trends)) warning("Cross-platform Trends CSV not found; do not infer directions")
if(!is.null(trends)) writecsv(trends,"04_Previous_Bulk_scRNA_Trends.csv")

wb <- openxlsx::createWorkbook()
add_sheet <- function(name,df) {
  openxlsx::addWorksheet(wb,name,gridLines=FALSE)
  hd <- openxlsx::createStyle(fgFill="#17324D",fontColour="#FFFFFF",
       textDecoration="bold",wrapText=TRUE)
  openxlsx::writeData(wb,name,df,withFilter=TRUE,headerStyle=hd)
  openxlsx::freezePane(wb,name,firstRow=TRUE)
  openxlsx::setColWidths(wb,name,cols=seq_along(df),widths=24)
  openxlsx::setRowHeights(wb,name,rows=1,heights=42)
}
readme <- data.frame(
  ITEM=c("Analysis status","CellChat centrality","Stage comparison limitation",
         "CellChat missing parameters","CIBERSORT provenance",
         "Cross-platform denominator","Cross-platform interpretation"),
  INTERPRETATION=c(
    "Read existing RDS and CSV only; NO re-fitting performed",
    "Weighted incoming/outgoing strength from @net$weight; ranks WITHIN stage only",
    "N of cell types and cells differs by stage: no between-stage inference test",
    "Original min.cells and DB release not verified; locate original analysis script",
    "CIBERSORT signature matrix, permutations, and mode require original run records",
    "Previously saved bulk/scRNA comparison normalizes four major immune lineages; not total tissue fractions",
    "Different cohorts/platforms; report concordant AND discordant directions descriptively"
  )
)
add_sheet("README",readme)
add_sheet("CellChat_Parameters",parameters)
add_sheet("CellChat_Strength",centrality)
add_sheet("CellChat_Networks",network)
if(!is.null(trends)) add_sheet("Prior_Bulk_scRNA_Trends",trends)
if(!is.null(samples)) add_sheet("Prior_Bulk_scRNA_Samples",samples)
file <- file.path(OUT,"CellChat_CIBERSORT_PreviousResults_Audit.xlsx")
openxlsx::saveWorkbook(wb,file,overwrite=TRUE)
cat("\nSAVED: ",file,"\n",sep="")
cat("Important: within-stage outgoing/incoming weight rank is NOT patient-level centrality significance.\n")

}

# ---- legacy_05_QC_scDblFinder_DecontX_Complete : Original independent post-hoc QC analysis ----
legacy_05_QC_scDblFinder_DecontX_Complete <- function() {
# ============================================================================
# PLOS ONE gastric scRNA-seq revision | Four-module independent QC audit
# 2026-10-21 | SOURCE OBJECTS ARE READ-ONLY
# 1) scDblFinder per original tissue sample (full cell atlas only)
# 2) DecontX ambient RNA evaluation (independent counts; no overwrite)
# 3) Five myeloid-state retention/marker/sample-support sensitivity audit
# 4) CellChat, Slingshot and CIBERSORT saved-parameter extraction
# Usage: in R with the Fig1_Original81779_Revised directory as working dir:
#        source(file.choose(), encoding = "UTF-8")
# Re-run the same script to resume completed samples.
# ============================================================================

# ------- SETTINGS: change only these if necessary ---------------------------
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
OUT <- file.path(ROOT, "QC_Methods_Audit_20261021", "FourModules_v1")
FULL_ATLAS_RDS <- ""  # Optional: absolute path to the exact 81,779-cell Seurat RDS
AUTO_INSTALL_AMBIENT <- TRUE  # Attempt Bioconductor decontX/celda if missing
RUN_DOUBLET <- TRUE
RUN_AMBIENT <- TRUE
RUN_SENSITIVITY <- TRUE
RUN_RECLUSTER_SENSITIVITY <- TRUE  # Sidecar clustering only; no original objects changed
RUN_METHOD_PARAMETERS <- TRUE
SAVE_AMBIENT_CORRECTED_COUNTS <- TRUE   # Separate, per-sample RDS: may use GB of disk
MIN_CELLS_PER_CAPTURE <- 80L
SEED <- 20261021L
# Descriptive sensitivity thresholds; DO NOT remove cells from published objects:
AMBIENT_SENSITIVITY_THRESHOLDS <- c(0.10, 0.20, 0.30)

if (basename(ROOT) != "Fig1_Original81779_Revised") {
  stop("Set working directory to .../Fig1_Original81779_Revised first. Current: ", ROOT)
}
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (p in c("01_Doublet", "02_Ambient", "03_Stability", "04_Methods")) {
  dir.create(file.path(OUT, p), recursive = TRUE, showWarnings = FALSE)
}
options(stringsAsFactors = FALSE, timeout = max(600, getOption("timeout", 60)))
set.seed(SEED)

logfile <- file.path(OUT, "00_Run_Log.txt")
logmsg <- function(...) {
  s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(s, "\n")
  cat(s, "\n", file = logfile, append = TRUE)
}
safe_write <- function(df, path) utils::write.csv(df, path, row.names = FALSE, na = "")
safe_read <- function(path) tryCatch(readRDS(path), error = function(e) NULL)
safe_rate <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_med <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
sanitize <- function(x) gsub("[^[:alnum:]_.-]", "_", as.character(x))
run_step <- function(name, fun) {
  logmsg("START ", name)
  tryCatch({
    z <- fun()
    logmsg("FINISH ", name)
    z
  }, error = function(e) {
    logmsg("ERROR ", name, ": ", conditionMessage(e))
    NULL
  })
}

required <- c("SeuratObject", "Matrix", "SingleCellExperiment", "SummarizedExperiment")
notfound <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(notfound)) stop("Missing essential packages: ", paste(notfound, collapse = ", "))

logmsg("BEGIN | directory=", ROOT)
logmsg("Output=", OUT)
logmsg("R=", R.version.string)

# ------- Load 3101 candidate cells (existing RDS, without modifications) ----
mac_path <- file.path(ROOT, "Fig4_Macrophage_Reanalysis", "Phase6_HighResolution",
                      "07_HighRes_Annotated_Candidates.rds")
if (!file.exists(mac_path)) stop("Cannot locate 3101-cell candidate RDS: ", mac_path)
if (exists("mac", envir = .GlobalEnv, inherits = FALSE) &&
    methods::is(get("mac", envir = .GlobalEnv), "Seurat") &&
    ncol(get("mac", envir = .GlobalEnv)) == 3101L) {
  macQC <- get("mac", envir = .GlobalEnv)
} else {
  macQC <- readRDS(mac_path)
}
stopifnot(ncol(macQC) == 3101L, "cluster_r1_0" %in% colnames(macQC@meta.data))
macCells <- colnames(macQC)
macMeta <- macQC@meta.data
macMeta$cell <- rownames(macMeta)
macMeta$cluster <- as.character(macMeta$cluster_r1_0)
macMeta$sample_audit <- if ("sample" %in% names(macMeta)) as.character(macMeta$sample) else
  if ("sample_uid" %in% names(macMeta)) as.character(macMeta$sample_uid) else NA_character_
macMeta$stage_audit <- if ("stage" %in% names(macMeta)) as.character(macMeta$stage) else NA_character_
logmsg("Macrophage candidate cells=", ncol(macQC), "; genes=", nrow(macQC))

# -------- Find the exact full 81779-cell Seurat object, never use 86422 ------
find_full_atlas <- function() {
  if (exists("Original81779_Fig1_Annotated_FINAL", envir = .GlobalEnv, inherits = FALSE)) {
    x <- get("Original81779_Fig1_Annotated_FINAL", envir = .GlobalEnv)
    if (methods::is(x, "Seurat") && ncol(x) == 81779L) {
      logmsg("Using existing in-memory Original81779_Fig1_Annotated_FINAL")
      return(x)
    }
  }
  roots <- unique(c(ROOT, dirname(ROOT), dirname(dirname(ROOT))))
  found <- unique(unlist(lapply(roots, function(x) {
    if (!dir.exists(x)) return(character())
    list.files(x, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE,
               ignore.case = TRUE)
  }), use.names = FALSE))
  exact <- if (nzchar(FULL_ATLAS_RDS)) FULL_ATLAS_RDS else character()
  wanted <- found[grepl("81779|Original.*FINAL|Fig1.*annotat|annotat.*Fig1", basename(found),
                        ignore.case = TRUE)]
  ordered <- unique(c(exact, wanted))
  safe_write(data.frame(path = found, bytes = file.info(found)$size,
                        preferred = found %in% ordered),
             file.path(OUT, "00_Full_Atlas_RDS_Inventory.csv"))
  if (!length(ordered)) {
    logmsg("No candidate full-atlas RDS found; set FULL_ATLAS_RDS in SETTINGS or load object into R.")
    return(NULL)
  }
  for (f in ordered) {
    if (!file.exists(f)) next
    logmsg("Inspecting candidate atlas RDS: ", f)
    obj <- tryCatch(readRDS(f), error = function(e) {
      logmsg("Cannot load: ", conditionMessage(e)); NULL
    })
    if (!is.null(obj) && methods::is(obj, "Seurat") && ncol(obj) == 81779L) {
      logmsg("Confirmed full atlas: ", f)
      return(obj)
    }
    logmsg("Not the 81779-cell Seurat atlas: ", f)
    rm(obj); gc(verbose = FALSE)
  }
  NULL
}

fullQC <- NULL
need_full <- RUN_DOUBLET || RUN_AMBIENT
if (need_full) fullQC <- run_step("find-full-atlas", find_full_atlas)
if (!is.null(fullQC)) {
  stopifnot(ncol(fullQC) == 81779L, "sample" %in% colnames(fullQC@meta.data))
  allSamples <- as.character(fullQC$sample)
  if (anyNA(allSamples) || any(!nzchar(allSamples))) {
    logmsg("Full-atlas sample has missing values; cannot perform per-sample QC safely.")
    fullQC <- NULL
  }
}

# ----------------------- Module 1: Doublet -------------------------------
run_doublet_module <- function() {
  dest <- file.path(OUT, "01_Doublet")
  if (is.null(fullQC)) {
    logmsg("SKIPPED formal doublet: exact full-atlas object unavailable")
    return(NULL)
  }
  if (!requireNamespace("scDblFinder", quietly = TRUE) ||
      !requireNamespace("BiocParallel", quietly = TRUE)) {
    logmsg("SKIPPED doublet: scDblFinder or BiocParallel not installed")
    return(NULL)
  }
  cn <- SeuratObject::GetAssayData(fullQC, assay = "RNA", layer = "counts")
  stopifnot(ncol(cn) == 81779L, identical(colnames(cn), colnames(fullQC)))
  indices <- split(seq_along(allSamples), allSamples)
  samples <- names(indices)
  safe_write(data.frame(sample = samples, total_cells = lengths(indices)),
             file.path(dest, "00_Sample_Input_Counts.csv"))
  logmsg("Formal scDblFinder starting; per-sample runs=", length(samples))
  for (i in seq_along(samples)) {
    sid <- samples[i]
    ix <- indices[[i]]
    fname <- file.path(dest, paste0("doublet_", sanitize(sid), ".rds"))
    if (file.exists(fname)) {
      logmsg("[", i, "/", length(samples), "] cached: ", sid)
      next
    }
    if (length(ix) < MIN_CELLS_PER_CAPTURE) {
      logmsg("[", i, "] SKIP too few cells: ", sid, " n=", length(ix))
      next
    }
    logmsg("[", i, "/", length(samples), "] scDblFinder: ", sid, " n=", length(ix))
    ok <- tryCatch({
      subcn <- cn[, ix, drop = FALSE]
      sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = subcn))
      set.seed(SEED + i)
      sce <- scDblFinder::scDblFinder(
        sce, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
      )
      cd <- as.data.frame(SummarizedExperiment::colData(sce))
      if (!all(c("scDblFinder.score", "scDblFinder.class") %in% names(cd))) {
        stop("scDblFinder missing expected results columns")
      }
      tab <- data.frame(cell = colnames(sce), sample = sid,
                        doublet_score = as.numeric(cd$scDblFinder.score),
                        doublet_class = as.character(cd$scDblFinder.class))
      saveRDS(tab, fname, compress = "gzip")
      logmsg("Saved ", sid, " | doubles=", sum(tab$doublet_class == "doublet", na.rm = TRUE))
      rm(sce, subcn, tab); gc(verbose = FALSE)
      TRUE
    }, error = function(e) {
      logmsg("FAILED sample ", sid, ": ", conditionMessage(e)); FALSE
    })
  }
  parts <- list.files(dest, pattern = "^doublet_.*\\.rds$", full.names = TRUE)
  if (!length(parts)) return(NULL)
  dt <- do.call(rbind, lapply(parts, readRDS))
  if (anyDuplicated(dt$cell)) stop("Doublet results have duplicated cell identifiers")
  safe_write(dt, file.path(dest, "01_All_Sample_Doublet_Calls.csv"))
  summary <- do.call(rbind, lapply(split(dt, dt$sample), function(x) {
    data.frame(sample = x$sample[1], tested = nrow(x),
               n_doublets = sum(x$doublet_class == "doublet", na.rm = TRUE),
               doublet_fraction = mean(x$doublet_class == "doublet", na.rm = TRUE),
               median_score = safe_med(x$doublet_score))
  }))
  safe_write(summary, file.path(dest, "02_Sample_Doublet_Summary.csv"))
  logmsg("scDblFinder complete files=", length(parts), " cells tested=", nrow(dt))
  dt
}
doublets <- NULL
if (RUN_DOUBLET) doublets <- run_step("module-1-doublet", run_doublet_module)

# ----------------------- Module 2: Ambient -------------------------------
ambient_engine <- function() {
  if (requireNamespace("decontX", quietly = TRUE)) return("decontX")
  if (requireNamespace("celda", quietly = TRUE)) return("celda")
  if (!AUTO_INSTALL_AMBIENT) return(NA_character_)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    tryCatch(utils::install.packages("BiocManager", repos = "https://cloud.r-project.org"),
             error = function(e) logmsg("BiocManager install error: ", conditionMessage(e)))
  }
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    for (p in c("decontX", "celda")) {
      logmsg("Attempting Bioconductor package installation: ", p)
      tryCatch(BiocManager::install(p, ask = FALSE, update = FALSE, quiet = TRUE),
               error = function(e) logmsg("Package install error for ", p, ": ", conditionMessage(e)))
      if (requireNamespace(p, quietly = TRUE)) return(p)
    }
  }
  NA_character_
}
run_ambient_module <- function() {
  dest <- file.path(OUT, "02_Ambient")
  # Record raw droplet sources if visible; matched droplet/background pairing
  # is NOT assumed just because a path contains 'raw'.
  audit_root <- sub("/GSE150290_fullcount_integration/.*$", "", ROOT)
  scan_dirs <- unique(c(if (dir.exists(audit_root)) audit_root else ROOT))
  candidates <- unique(unlist(lapply(scan_dirs, function(d) {
    if (!dir.exists(d)) return(character())
    list.files(d, pattern = "raw_feature_bc_matrix|raw_matrix|unfiltered",
               recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))
  safe_write(data.frame(possible_raw_droplet_source = candidates),
             file.path(dest, "00_Raw_Droplet_File_Check.csv"))
  if (is.null(fullQC)) {
    logmsg("SKIPPED formal ambient: exact full atlas unavailable")
    return(NULL)
  }
  eng <- ambient_engine()
  writeLines(if (is.na(eng)) "NOT_INSTALLED" else eng,
             file.path(dest, "00_Ambient_Engine.txt"))
  if (is.na(eng)) {
    logmsg("SKIPPED formal ambient: decontX/celda package unavailable; see package log")
    return(NULL)
  }
  logmsg("Ambient engine=", eng, "; using full filtered sample counts; no raw droplet correction claimed")
  cn <- SeuratObject::GetAssayData(fullQC, assay = "RNA", layer = "counts")
  indices <- split(seq_along(allSamples), allSamples)
  marker_genes <- c("GAST", "SST", "CHGA", "CHGB", "EPCAM", "KRT8",
                    "KRT18", "KRT19", "PIGR", "PHGR1", "CD3D", "CD3E",
                    "LST1", "TYROBP", "FCER1G", "C1QA", "SPP1", "TREM2")
  for (i in seq_along(indices)) {
    sid <- names(indices)[i]
    ix <- indices[[i]]
    fname <- file.path(dest, paste0("ambient_", sanitize(sid), ".rds"))
    if (file.exists(fname)) {logmsg("Ambient cached: ", sid); next}
    if (length(ix) < MIN_CELLS_PER_CAPTURE) {logmsg("Ambient skip small sample: ", sid); next}
    logmsg("[", i, "/", length(indices), "] DecontX: ", sid, " n=", length(ix))
    tryCatch({
      raw <- cn[, ix, drop = FALSE]
      sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = raw))
      set.seed(SEED + 1000L + i)
      # Per sample, z=NULL: learn expression clusters from that sample.
      if (eng == "decontX") {
        res <- decontX::decontX(sce, maxIter = 200, seed = SEED + 1000L + i,
                                verbose = FALSE)
        cleaned <- decontX::decontXcounts(res)
      } else {
        res <- celda::decontX(sce, maxIter = 200, seed = SEED + 1000L + i,
                              verbose = FALSE)
        cleaned <- celda::decontXcounts(res)
      }
      contam <- as.numeric(SummarizedExperiment::colData(res)$decontX_contamination)
      if (length(contam) != ncol(raw)) stop("Unexpected contamination vector length")
      genes <- intersect(marker_genes, rownames(raw))
      before <- as.matrix(raw[genes, , drop = FALSE])
      after <- as.matrix(cleaned[genes, , drop = FALSE])
      marker_tab <- data.frame(
        cell = rep(colnames(raw), each = length(genes)),
        sample = sid,
        gene = rep(genes, times = ncol(raw)),
        raw_count = as.numeric(before),
        decont_count = as.numeric(after)
      )
      metrics <- data.frame(cell = colnames(raw), sample = sid,
                            contamination = contam)
      saveRDS(list(metrics = metrics, markers = marker_tab, engine = eng,
                   corrected_counts = if (SAVE_AMBIENT_CORRECTED_COUNTS) cleaned else NULL),
              fname, compress = "gzip")
      logmsg("Saved ambient ", sid, " | median contamination=", round(safe_med(contam), 4))
      rm(raw, sce, res, cleaned, before, after, marker_tab, metrics); gc(verbose = FALSE)
    }, error = function(e) logmsg("FAILED ambient sample ", sid, ": ", conditionMessage(e)))
  }
  pieces <- list.files(dest, pattern = "^ambient_.*\\.rds$", full.names = TRUE)
  if (!length(pieces)) return(NULL)
  # Save and merge only summaries/markers; do not merge huge corrected matrices.
  mets <- list(); marks <- list(); manifests <- list()
  for (f in pieces) {
    a <- readRDS(f)
    mets[[length(mets) + 1L]] <- a$metrics
    marks[[length(marks) + 1L]] <- a$markers
    manifests[[length(manifests) + 1L]] <- data.frame(
      file = basename(f), sample = a$metrics$sample[1], n_cells = nrow(a$metrics),
      engine = a$engine, corrected_counts_saved = !is.null(a$corrected_counts)
    )
    rm(a); gc(verbose = FALSE)
  }
  am <- do.call(rbind, mets)
  mr <- do.call(rbind, marks)
  safe_write(am, file.path(dest, "01_All_Cell_Contamination.csv"))
  safe_write(mr, file.path(dest, "02_Marker_Before_After.csv"))
  safe_write(do.call(rbind, manifests), file.path(dest, "03_Corrected_Counts_Manifest.csv"))
  logmsg("Ambient results cells=", nrow(am), "; samples=", length(pieces))
  list(metrics = am, markers = mr)
}
ambient <- NULL
if (RUN_AMBIENT) ambient <- run_step("module-2-ambient", run_ambient_module)

# --------- Module 3: sensitivity / retention audit, not reclustering ---------
run_stability_module <- function() {
  dest <- file.path(OUT, "03_Stability")
  selected <- c("0", "1", "3", "7", "11")
  states <- c("0" = "FCN1+ monocyte-like", "1" = "SPP1+ macrophages",
              "3" = "C1Q+ macrophages", "7" = "APOE+/TREM2+ macrophages",
              "11" = "IL1B+ inflammatory myeloid")
  d <- macMeta[, c("cell", "cluster", "sample_audit", "stage_audit"), drop = FALSE]
  d$state <- unname(states[d$cluster])
  d$selected <- d$cluster %in% selected
  d$doublet_class <- NA_character_
  d$doublet_score <- NA_real_
  d$contamination <- NA_real_
  if (!is.null(doublets)) {
    mm <- match(d$cell, doublets$cell)
    d$doublet_class <- doublets$doublet_class[mm]
    d$doublet_score <- doublets$doublet_score[mm]
  }
  if (!is.null(ambient)) {
    mm <- match(d$cell, ambient$metrics$cell)
    d$contamination <- ambient$metrics$contamination[mm]
  }
  # A missing result means UNTESTED, NOT singlet/low contamination.
  safe_write(d, file.path(dest, "01_All3101_Audit_Mapping.csv"))
  mapping_summary <- data.frame(
    total_candidates = nrow(d),
    doublet_tested = sum(!is.na(d$doublet_class)),
    ambient_tested = sum(!is.na(d$contamination)),
    selected_cells = sum(d$selected),
    selected_doublet_tested = sum(d$selected & !is.na(d$doublet_class)),
    selected_ambient_tested = sum(d$selected & !is.na(d$contamination))
  )
  safe_write(mapping_summary, file.path(dest, "00_Data_Coverage.csv"))
  logmsg("QC mapping: ", paste(names(mapping_summary), mapping_summary[1, ], collapse = "; "))
  # marker identity sensitivity before/after doublet calls; no label is changed.
  cn <- SeuratObject::GetAssayData(macQC, assay = "RNA", layer = "counts")
  marker_sets <- list(
    myeloid = c("LST1", "TYROBP", "FCER1G", "AIF1", "CD68", "CSF1R"),
    epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "PIGR", "PHGR1"),
    tcell = c("CD3D", "CD3E", "CD2", "TRBC1"),
    neuroendocrine = c("GAST", "SST", "CHGA", "CHGB")
  )
  for (nm in names(marker_sets)) {
    genes <- intersect(marker_sets[[nm]], rownames(cn))
    d[[paste0(nm, "_two_markers")]] <- if (length(genes) < 2L) {
      rep(NA, nrow(d))
    } else as.vector(Matrix::colSums(cn[genes, , drop = FALSE] > 0) >= 2)
  }
  d$mixed_myeloid_epithelial <- with(d, myeloid_two_markers & epithelial_two_markers)
  d$mixed_myeloid_tcell <- with(d, myeloid_two_markers & tcell_two_markers)
  d$predicted_doublet <- ifelse(is.na(d$doublet_class), NA,
                                d$doublet_class == "doublet")
  scenarios <- list(baseline = rep(TRUE, nrow(d)))
  # Only compare an alternative when all selected cells have the required QC.
  if (all(!is.na(d$doublet_class[d$selected]))) {
    scenarios$exclude_predicted_doublets <- d$doublet_class != "doublet"
  } else logmsg("NO complete doublet sensitivity scenario: selected cells partially untested")
  if (all(!is.na(d$contamination[d$selected]))) {
    for (threshold in AMBIENT_SENSITIVITY_THRESHOLDS) {
      scenarios[[paste0("ambient_fraction_le_", threshold)]] <-
        d$contamination <= threshold
      if ("exclude_predicted_doublets" %in% names(scenarios)) {
        scenarios[[paste0("exclude_doublets_ambient_le_", threshold)]] <-
          (d$doublet_class != "doublet") & (d$contamination <= threshold)
      }
    }
  } else logmsg("NO complete ambient sensitivity scenario: selected cells partially untested")
  stats <- list(); sample_stats <- list()
  for (nm in names(scenarios)) {
    valid <- scenarios[[nm]]
    valid[is.na(valid)] <- FALSE
    for (cl in selected) {
      idx <- which(d$cluster == cl)
      kept <- idx[valid[idx]]
      st <- d$sample_audit[kept]
      st <- st[!is.na(st) & nzchar(st)]
      tb <- table(st)
      stats[[length(stats) + 1L]] <- data.frame(
        scenario = nm, cluster = cl, state = states[[cl]],
        original_cells = length(idx), retained_cells = length(kept),
        retention_rate = length(kept) / length(idx),
        n_tissue_samples = length(tb), n_samples_at_least5_cells = sum(tb >= 5),
        myeloid_flag_rate = safe_rate(d$myeloid_two_markers[kept]),
        epithelial_flag_rate = safe_rate(d$epithelial_two_markers[kept]),
        myeloid_epithelial_mixed_rate = safe_rate(d$mixed_myeloid_epithelial[kept]),
        tcell_flag_rate = safe_rate(d$tcell_two_markers[kept]),
        median_doublet_score = safe_med(d$doublet_score[kept]),
        median_contamination = safe_med(d$contamination[kept])
      )
      if (length(tb)) sample_stats[[length(sample_stats) + 1L]] <- data.frame(
        scenario = nm, cluster = cl, sample = names(tb), n_cells = as.integer(tb))
    }
  }
  safe_write(do.call(rbind, stats), file.path(dest, "02_State_Retention_Sensitivity.csv"))
  if (length(sample_stats)) safe_write(do.call(rbind, sample_stats),
                                      file.path(dest, "03_Tissue_Sample_Support.csv"))
  # Expression-marker retention: what remains in each candidate state.
  genes <- c("FCN1", "VCAN", "S100A8", "IL1B", "CCL3", "CCL4",
             "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "APOC1",
             "SPP1", "GPNMB", "PLA2G7", "GAST", "SST", "EPCAM", "KRT19", "CD3D")
  genes <- intersect(genes, rownames(cn))
  gene_summary <- list()
  for (nm in names(scenarios)) for (cl in selected) {
    keep <- which(d$cluster == cl & (scenarios[[nm]] %in% TRUE))
    if (!length(keep)) next
    v <- cn[genes, keep, drop = FALSE]
    gene_summary[[length(gene_summary) + 1L]] <- data.frame(
      scenario = nm, cluster = cl, state = states[[cl]], gene = genes,
      detected_fraction = as.numeric(Matrix::rowMeans(v > 0)),
      mean_raw_UMI = as.numeric(Matrix::rowMeans(v))
    )
  }
  if (length(gene_summary)) safe_write(do.call(rbind, gene_summary),
                                      file.path(dest, "04_State_Marker_Retention.csv"))
  if (!is.null(ambient)) {
    subset_markers <- ambient$markers[ambient$markers$cell %in% d$cell[d$selected], , drop = FALSE]
    if (nrow(subset_markers)) {
      safe_write(subset_markers, file.path(dest, "05_Selected_State_Ambient_Markers.csv"))
      # Avoid zero division; mean change only for cells with positive original counts.
      before_after <- aggregate(cbind(raw_count, decont_count) ~ gene, subset_markers, sum)
      before_after$fraction_removed <- ifelse(before_after$raw_count > 0,
                                               1 - before_after$decont_count / before_after$raw_count,
                                               NA_real_)
      safe_write(before_after, file.path(dest, "06_Overall_Marker_Before_After.csv"))
      state_lookup <- d[, c("cell", "state"), drop = FALSE]
      ss <- merge(subset_markers, state_lookup, by = "cell", all.x = TRUE)
      ss <- ss[!is.na(ss$state), , drop = FALSE]
      if (nrow(ss)) {
        per_state <- aggregate(cbind(raw_count, decont_count) ~ state + gene, ss, sum)
        per_state$fraction_removed <- ifelse(per_state$raw_count > 0,
          1 - per_state$decont_count / per_state$raw_count, NA_real_)
        safe_write(per_state, file.path(dest, "07_State_Marker_Before_After.csv"))
      }
    }
  }
  # Optional sidecar reclustering at resolution 1.0. This checks whether
  # retaining high-confidence cells preserves the original partition structure;
  # it does NOT establish lineage identity or biological reproducibility.
  if (RUN_RECLUSTER_SENSITIVITY && requireNamespace("Seurat", quietly = TRUE)) {
    candidate_scenarios <- intersect(
      c("exclude_predicted_doublets", "exclude_doublets_ambient_le_0.2"),
      names(scenarios)
    )
    if (length(candidate_scenarios) &&
        all(!is.na(d$doublet_class)) &&
        (is.null(ambient) || all(!is.na(d$contamination)))) {
      for (nm in candidate_scenarios) {
        fn <- file.path(dest, paste0("05_Recluster_", sanitize(nm), ".rds"))
        if (file.exists(fn)) {logmsg("Reclustering cached: ", nm); next}
        keep <- which(scenarios[[nm]] %in% TRUE)
        if (length(keep) < 300) {logmsg("Reclustering skipped; too few cells: ", nm); next}
        logmsg("SIDE-CAR reclustering: ", nm, " n=", length(keep))
        tryCatch({
          subobj <- subset(macQC, cells = d$cell[keep])
          subobj <- SeuratObject::`DefaultAssay<-`(subobj, value = "RNA")
          subobj <- Seurat::NormalizeData(subobj, verbose = FALSE)
          subobj <- Seurat::FindVariableFeatures(subobj, nfeatures = 2500, verbose = FALSE)
          hv <- SeuratObject::VariableFeatures(subobj)
          hv <- hv[!grepl("^(MT-|RPL|RPS)", hv)]
          if (length(hv) < 200) stop("Not enough variable features after filtering")
          subobj <- SeuratObject::`VariableFeatures<-`(subobj, value = hv)
          subobj <- Seurat::ScaleData(subobj, features = hv, verbose = FALSE)
          subobj <- Seurat::RunPCA(subobj, features = hv, npcs = 20, verbose = FALSE)
          subobj <- Seurat::FindNeighbors(subobj, reduction = "pca", dims = 1:20,
                                       verbose = FALSE)
          subobj <- Seurat::FindClusters(subobj, resolution = 1.0, verbose = FALSE)
          new <- as.character(SeuratObject::Idents(subobj))
          old <- d$cluster[match(colnames(subobj), d$cell)]
          ct <- as.data.frame(table(original_cluster = old, new_cluster = new))
          ct <- ct[ct$Freq > 0, , drop = FALSE]
          safe_write(ct, sub("\\.rds$", "_CrossTab.csv", fn))
          ari <- if (requireNamespace("mclust", quietly = TRUE)) {
            mclust::adjustedRandIndex(old, new)
          } else NA_real_
          saveRDS(list(scenario = nm, retained_n = length(keep),
                       orig_clusters = length(unique(old)),
                       new_clusters = length(unique(new)),
                       adjusted_Rand_index = ari,
                       cell_cluster_assignments = data.frame(
                         cell = colnames(subobj), original_cluster = old, new_cluster = new)),
                  fn, compress = "gzip")
          logmsg("Sidecar reclustering finished: ", nm, " ARI=", ari)
          rm(subobj); gc(verbose = FALSE)
        }, error = function(e) logmsg("Sidecar reclustering FAILED: ", nm,
                                      " | ", conditionMessage(e)))
      }
    } else logmsg("Sidecar reclustering skipped: requires complete 3101-cell QC coverage")
  }
  logmsg("Sensitivity audit complete: retention and optional sidecar reclustering.")
  invisible(mapping_summary)
}
if (RUN_SENSITIVITY) run_step("module-3-state-sensitivity", run_stability_module)

# ----------------------- Module 4: saved parameters -----------------------
run_parameter_module <- function() {
  dest <- file.path(OUT, "04_Methods")
  pkg <- c("Seurat", "SeuratObject", "scDblFinder", "decontX", "celda",
           "CellChat", "slingshot", "edgeR", "limma", "fgsea")
  versions <- vapply(pkg, function(x) {
    if (requireNamespace(x, quietly = TRUE)) as.character(utils::packageVersion(x)) else NA_character_
  }, character(1))
  safe_write(data.frame(package = pkg, version = versions),
             file.path(dest, "00_Package_Versions.csv"))
  if ("commands" %in% methods::slotNames(macQC)) {
    con <- file(file.path(dest, "01_Macrophage_Seurat_Commands.txt"), "wt")
    for (nm in names(macQC@commands)) {
      cmd <- macQC@commands[[nm]]
      writeLines(paste0("=== ", nm, " ==="), con)
      params <- tryCatch(methods::slot(cmd, "params"), error = function(e) NULL)
      writeLines(capture.output(str(params, max.level = 3)), con)
    }
    close(con)
  }
  ccdir <- file.path(ROOT, "Fig2_CellChat_Revised")
  ccfile <- file.path(ccdir, paste0("CellChat_", c("NAG", "CAG", "IM", "GC"), "_primary.rds"))
  chatstats <- list()
  con <- file(file.path(dest, "02_CellChat_Saved_Options.txt"), "wt")
  for (f in ccfile) {
    stage <- sub("^CellChat_(.*)_primary\\.rds$", "\\1", basename(f))
    if (!file.exists(f)) {writeLines(paste(stage, "NOT FOUND"), con); next}
    cc <- tryCatch(readRDS(f), error = function(e) {writeLines(conditionMessage(e), con); NULL})
    if (is.null(cc)) next
    writeLines(paste("===", stage, "==="), con)
    op <- tryCatch(cc@options, error = function(e) NULL)
    writeLines(capture.output(str(op, max.level = 4)), con)
    m <- tryCatch(cc@net$count, error = function(e) NULL)
    weight <- tryCatch(cc@net$weight, error = function(e) NULL)
    chatstats[[length(chatstats) + 1L]] <- data.frame(
      stage = stage, n_cells = nrow(cc@meta),
      n_celltypes = if (!is.null(m)) nrow(m) else NA_integer_,
      inferred_interaction_count = if (!is.null(m)) sum(m, na.rm = TRUE) else NA_real_,
      inferred_weight_sum = if (!is.null(weight)) sum(weight, na.rm = TRUE) else NA_real_,
      database_name = if (!is.null(op$database)) as.character(op$database)[1] else NA_character_
    )
    rm(cc); gc(verbose = FALSE)
  }
  close(con)
  if (length(chatstats)) safe_write(do.call(rbind, chatstats),
                                   file.path(dest, "03_CellChat_Network_Summary.csv"))
  # Read existing Slingshot inputs / files ONLY. Never infer root from a figure.
  audit_root <- sub("/GSE150290_fullcount_integration/.*$", "", ROOT)
  candidate_dirs <- c(file.path(ROOT, "Fig4_Macrophage_Reanalysis"),
                      file.path(audit_root, "Fig4_Final"))
  sf <- unique(unlist(lapply(candidate_dirs, function(d) {
    if (!dir.exists(d)) return(character())
    list.files(d, pattern = "slingshot|pseudotime|trajectory|lineage|Fig4C_Main",
               full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))
  safe_write(data.frame(path = sf, file_size = file.info(sf)$size),
             file.path(dest, "04_Slingshot_File_Inventory.csv"))
  sc_rds <- sf[grepl("\\.rds$", sf, ignore.case = TRUE)]
  con <- file(file.path(dest, "05_Slingshot_Object_Inspection.txt"), "wt")
  if (!length(sc_rds)) {
    writeLines("NO saved Slingshot RDS identifiable by filename. Root and fitting parameters NOT VERIFIED.", con)
  } else for (f in sc_rds) {
    writeLines(paste("===", f, "==="), con)
    size <- file.info(f)$size
    if (is.na(size) || size > 400 * 1024^2) {
      writeLines("SKIPPED (>400MB; object not auto-loaded)", con); next
    }
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x)) {writeLines("Could not load", con); next}
    writeLines(capture.output({
      print(class(x)); print(names(x));
      print(tryCatch(x@metadata, error = function(e) "No S4 metadata"))
      if (requireNamespace("slingshot", quietly = TRUE)) {
        print(tryCatch(length(slingshot::slingLineages(x)), error = function(e) NA_integer_))
      }
    }), con)
    rm(x); gc(verbose = FALSE)
  }
  close(con)
  # CIBERSORT data and potential configuration: do not fabricate LM22/permutations.
  d3 <- file.path(ROOT, "Fig3_ABC_Revised")
  cb <- if (dir.exists(d3)) list.files(d3, recursive = TRUE, full.names = TRUE) else character()
  cb <- cb[grepl("cibersort|LM22|signature|bulk|Fig3D_|Fig3E_", basename(cb), ignore.case = TRUE)]
  safe_write(data.frame(path = cb, bytes = file.info(cb)$size),
             file.path(dest, "06_CIBERSORT_File_Inventory.csv"))
  con <- file(file.path(dest, "07_CIBERSORT_File_Headers.txt"), "wt")
  for (f in cb) {
    if (!grepl("\\.(csv|tsv|txt|R|Rmd)$", f, ignore.case = TRUE)) next
    writeLines(paste("===", f, "==="), con)
    z <- tryCatch(readLines(f, n = 12L, warn = FALSE), error = function(e) "READ FAILED")
    writeLines(z, con)
  }
  close(con)
  notes <- c(
    "CellChat: stored stage-specific options are audited. Do not infer original min.cells or population.size from the later sensitivity script.",
    "Slingshot: unless start.clus/root parameters are found in an actual saved fitting object/script, mark them UNVERIFIED.",
    "CIBERSORT: unless mixture/signature/permutation run file exists, LM22, absolute/relative mode and permutation count are UNVERIFIED.",
    "Ambient: per-sample DecontX without raw droplets estimates background internally; do not describe this as SoupX/raw-empty-droplet estimation.",
    "Doublet: scDblFinder doublet calls use the retained full-atlas cells per capture; excluded low-quality droplets are not retrospectively represented.",
    "Stability: cell retention after computational filtering is sensitivity, not proof of independent reclustering or cell-lineage identity."
  )
  writeLines(notes, file.path(dest, "08_Methods_Caveats.txt"))
  logmsg("Method parameter inventory complete")
  invisible(TRUE)
}
if (RUN_METHOD_PARAMETERS) run_step("module-4-parameters", run_parameter_module)

# -------------------------------- SUMMARY ----------------------------------
eligible_samples <- if (!is.null(fullQC)) sum(table(allSamples) >= MIN_CELLS_PER_CAPTURE) else NA_integer_
n_doublet_samples <- length(list.files(file.path(OUT, "01_Doublet"),
                                     pattern = "^doublet_.*\\.rds$"))
n_ambient_samples <- length(list.files(file.path(OUT, "02_Ambient"),
                                     pattern = "^ambient_.*\\.rds$"))
classify_completion <- function(n) {
  if (is.na(eligible_samples)) return("NO_FULL_ATLAS")
  if (n == 0L) return("NOT_RUN")
  if (n < eligible_samples) return("PARTIAL")
  "COMPLETE"
}
status <- data.frame(
  module = c("formal_scDblFinder", "formal_DecontX", "five_state_sensitivity", "methods_parameters"),
  output_directory = file.path(OUT, c("01_Doublet", "02_Ambient", "03_Stability", "04_Methods")),
  status = c(classify_completion(n_doublet_samples),
             classify_completion(n_ambient_samples),
             if (file.exists(file.path(OUT, "03_Stability", "02_State_Retention_Sensitivity.csv")))
               "SUMMARY_WRITTEN_CHECK_COVERAGE" else "NOT_RUN",
             if (file.exists(file.path(OUT, "04_Methods", "00_Package_Versions.csv")))
               "INVENTORY_WRITTEN" else "NOT_RUN"),
  tested_samples = c(n_doublet_samples, n_ambient_samples, NA_integer_, NA_integer_),
  expected_samples = c(eligible_samples, eligible_samples, NA_integer_, NA_integer_)
)
safe_write(status, file.path(OUT, "99_Module_Status.csv"))
logmsg("ALL MODULES FINISHED OR RECORDED AS SKIPPED. ",
       "See 99_Module_Status.csv and 00_Run_Log.txt. No original data was overwritten.")
print(status)
cat("\nResults directory: ", OUT, "\n", sep = "")

}

# ---- legacy_06_Fig4_FourChecks : QC 4 checks ----
legacy_06_Fig4_FourChecks <- function() {
# ============================================================================
# PLOS ONE Fig4 | Four complementary robustness checks | 2026-10-21
# 1. State-marker enrichment and expression sensitivity
# 2. APOE+/TREM2+ identity, mixed-lineage and cross-sample support
# 3. Base-R adjusted Rand index (ARI), original vs QC-reclustered partitions
# 4. Fig4 functional GO and Slingshot sensitivity (with explicit limitations)
# Reads outputs of PLOS_GC_QC_FourModules_20261021.R; does NOT overwrite source.
# Run from the Fig1_Original81779_Revised folder in R: source(file.choose())
# ============================================================================

ROOT <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
stopifnot(basename(ROOT) == 'Fig1_Original81779_Revised')
QC <- file.path(ROOT, 'QC_Methods_Audit_20261021', 'FourModules_v1')
OUT <- file.path(QC, '05_Fig4_Validation_v1')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
STOP_ON_MISSING_QC <- TRUE
RUN_DE_AND_GO <- TRUE             # Rerun five-state DEG and GO-BP after QC filtering
RUN_SLINGSHOT_REFIT <- TRUE      # Secondary Slingshot fitting with saved embedding
SEED <- 20261021L
set.seed(SEED)

logfile <- file.path(OUT, '00_Fig4_Validation_Log.txt')
logmsg <- function(...) {
  x <- paste(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), paste0(..., collapse = ''), sep = ' | ')
  cat(x, '\n')
  cat(x, '\n', file = logfile, append = TRUE)
}
writecsv <- function(x, fn) utils::write.csv(x, file.path(OUT, fn), row.names = FALSE, na = '')
read_checked <- function(fn, required = TRUE) {
  if (!file.exists(fn)) {
    if (required && STOP_ON_MISSING_QC) stop('Missing required QC file: ', fn)
    logmsg('MISSING ', fn)
    return(NULL)
  }
  utils::read.csv(fn, stringsAsFactors = FALSE, check.names = FALSE)
}
safe_med <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
safe_mean <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
choose2 <- function(x) x * (x - 1) / 2
safe_ratio <- function(a,b) ifelse(b > 0, a/b, NA_real_)
ari_calc <- function(x, y) {
  stopifnot(length(x) == length(y))
  good <- !(is.na(x) | is.na(y)); x <- as.character(x[good]); y <- as.character(y[good])
  if (length(x) < 3) return(NA_real_)
  tt <- table(x, y)
  total <- choose2(sum(tt)); if (total <= 0) return(NA_real_)
  a <- sum(choose2(tt)); b <- sum(choose2(rowSums(tt)))
  c <- sum(choose2(colSums(tt))); expected <- b*c/total
  denominator <- (b+c)/2 - expected
  if (abs(denominator) < 1e-12) return(if (abs(a - expected) < 1e-12) 1 else NA_real_)
  (a-expected)/denominator
}

states <- c('0'='FCN1+ monocyte-like', '1'='SPP1+ macrophages',
            '3'='C1Q+ macrophages', '7'='APOE+/TREM2+ macrophages',
            '11'='IL1B+ inflammatory myeloid')
key_genes <- list(
  '0'=c('FCN1','VCAN','S100A8'),
  '1'=c('SPP1','GPNMB','PLA2G7'),
  '3'=c('C1QA','C1QB','C1QC'),
  '7'=c('APOE','TREM2','APOC1'),
  '11'=c('IL1B','CCL3','CCL4'))
selected <- names(states)
scenario_names <- c('baseline','exclude_predicted_doublets','exclude_doublets_ambient_le_0.2')

logmsg('START | Fig4 robustness | root=', ROOT)
logmsg('All original files and figures are READ-ONLY')

mapping <- read_checked(file.path(QC,'03_Stability','01_All3101_Audit_Mapping.csv'))
retention <- read_checked(file.path(QC,'03_Stability','02_State_Retention_Sensitivity.csv'))
pre_summary <- read_checked(file.path(QC,'03_Stability','04_State_Marker_Retention.csv'))
stopifnot(!anyDuplicated(mapping$cell), nrow(mapping) == 3101L)
stopifnot(all(c('cell','cluster','doublet_class','contamination') %in% names(mapping)))
stopifnot(all(c('scenario','cluster','gene','detected_fraction','mean_raw_UMI') %in% names(pre_summary)))
mapping$cluster <- as.character(mapping$cluster)
pre_summary$cluster <- as.character(pre_summary$cluster)
retention$cluster <- as.character(retention$cluster)
stopifnot(!anyNA(mapping$cell), !anyNA(mapping$cluster))
if (anyNA(mapping$doublet_class) || anyNA(mapping$contamination))
  stop('Some QC calls are missing; do not silently treat untested cells as low-risk.')

scenarios <- list(
  baseline = setNames(rep(TRUE, nrow(mapping)), mapping$cell),
  exclude_predicted_doublets = setNames(mapping$doublet_class != 'doublet', mapping$cell),
  exclude_doublets_ambient_le_0.2 = setNames(
    (mapping$doublet_class != 'doublet') & (mapping$contamination <= 0.2), mapping$cell))

# ---- 1: marker persistence from per-state existing sensitivity summaries ----
logmsg('CHECK1: Marker persistence')
comparisons <- list()
for (cl in selected) for (gene in key_genes[[cl]]) {
  z <- pre_summary[pre_summary$cluster == cl & pre_summary$gene == gene &
                   pre_summary$scenario %in% scenario_names, , drop=FALSE]
  base <- z[z$scenario == 'baseline', , drop=FALSE]
  if (!nrow(base)) next
  for (sc in scenario_names) {
    zz <- z[z$scenario == sc, , drop=FALSE]
    if (!nrow(zz)) next
    comparisons[[length(comparisons)+1L]] <- data.frame(
      cluster=cl,state=states[[cl]],gene=gene,scenario=sc,
      baseline_detected_fraction=base$detected_fraction[1],
      current_detected_fraction=zz$detected_fraction[1],
      delta_detected_fraction=zz$detected_fraction[1]-base$detected_fraction[1],
      baseline_mean_raw_UMI=base$mean_raw_UMI[1],
      current_mean_raw_UMI=zz$mean_raw_UMI[1],
      mean_raw_UMI_ratio=safe_ratio(zz$mean_raw_UMI[1],base$mean_raw_UMI[1]))
  }
}
if (length(comparisons)) {
  c1 <- do.call(rbind,comparisons)
  writecsv(c1, '01_CoreMarker_BeforeAfter_Filter.csv')
  print(c1[c1$cluster=='7',,drop=FALSE],row.names=FALSE)
} else logmsg('WARNING: no core markers matched existing 04_State_Marker_Retention.csv')

# ---- 2: marker enrichment / identity (also after DecontX correction) -------
logmsg('CHECK2: APOE/TREM2 identity and mixed-lineage signals')
macfile <- file.path(ROOT,'Fig4_Macrophage_Reanalysis','Phase6_HighResolution',
                     '07_HighRes_Annotated_Candidates.rds')
stopifnot(file.exists(macfile))
if (!requireNamespace('SeuratObject',quietly=TRUE) ||
    !requireNamespace('Matrix',quietly=TRUE)) stop('SeuratObject/Matrix are required')
mac <- readRDS(macfile)
stopifnot(ncol(mac)==3101L, all(mapping$cell %in% colnames(mac)))
ct <- SeuratObject::GetAssayData(mac, assay='RNA', layer='counts')

identity_genes <- unique(c(unlist(key_genes),
  'LST1','TYROBP','FCER1G','AIF1','CD68','CSF1R',
  'EPCAM','KRT8','KRT18','KRT19','PIGR','PHGR1',
  'CD3D','CD3E','GAST','SST','CHGA','CHGB'))
identity_genes <- intersect(identity_genes,rownames(ct))
# Enrichment relative to other candidates, with exactly the same QC inclusion.
identity_tbl <- list()
for (sc in scenario_names) {
  keep <- names(scenarios[[sc]])[scenarios[[sc]]]
  keep <- intersect(keep, colnames(ct))
  for (cl in selected) {
    inside <- intersect(keep, mapping$cell[mapping$cluster == cl])
    outside <- setdiff(keep, inside)
    if (!length(inside) || !length(outside)) next
    a <- ct[identity_genes, inside, drop=FALSE]
    b <- ct[identity_genes, outside, drop=FALSE]
    ai <- as.numeric(Matrix::rowMeans(a > 0)); bi <- as.numeric(Matrix::rowMeans(b > 0))
    ma <- as.numeric(Matrix::rowMeans(a)); mb <- as.numeric(Matrix::rowMeans(b))
    identity_tbl[[length(identity_tbl)+1L]] <- data.frame(
      scenario=sc, cluster=cl,state=states[[cl]],gene=identity_genes,
      n_state_cells=length(inside),n_other_cells=length(outside),
      pct_state=ai,pct_other_candidates=bi,delta_detection=ai-bi,
      mean_raw_UMI_state=ma,mean_raw_UMI_other=mb,
      log2_mean_UMI_ratio=log2((ma+0.01)/(mb+0.01)))
  }
}
identity_tbl <- do.call(rbind, identity_tbl)
writecsv(identity_tbl,'02_Marker_Enrichment_vs_OtherCandidates.csv')
print(identity_tbl[identity_tbl$cluster == '7' & identity_tbl$gene %in% key_genes[['7']],
                   c('scenario','gene','n_state_cells','pct_state','pct_other_candidates',
                     'delta_detection','log2_mean_UMI_ratio')],row.names=FALSE)

# Independently examine per-sample support, including dominant sample.
sample_tbl <- list()
for (sc in scenario_names) for (cl in selected) {
  mask <- mapping$cluster == cl & as.logical(scenarios[[sc]][mapping$cell])
  z <- mapping[mask,,drop=FALSE]
  if (!nrow(z)) next
  sample_col <- if ('sample_audit' %in% names(z)) 'sample_audit' else if ('sample' %in% names(z)) 'sample' else NA_character_
  if (is.na(sample_col)) stop('No tissue sample identifier in QC mapping')
  counts <- sort(table(z[[sample_col]]),decreasing=TRUE)
  sample_tbl[[length(sample_tbl)+1L]] <- data.frame(
    scenario=sc,cluster=cl,state=states[[cl]],n_cells=nrow(z),
    n_tissue_samples=length(counts),n_samples_ge5=sum(counts >=5),
    dominant_sample=names(counts)[1],dominant_n=as.numeric(counts[1]),
    dominant_fraction=as.numeric(counts[1])/nrow(z))
}
sample_tbl <- do.call(rbind,sample_tbl)
writecsv(sample_tbl,'03_State_TissueSupport_Concentration.csv')
print(sample_tbl[sample_tbl$cluster %in% c('7','11'),,drop=FALSE],row.names=FALSE)

# Summarize mixed signals, flag != confirmed doublet or ambient contamination.
mixed_tbl <- retention[retention$scenario %in% scenario_names &
                       retention$cluster %in% selected,
  intersect(c('scenario','cluster','state','original_cells','retained_cells',
     'retention_rate','myeloid_flag_rate','epithelial_flag_rate',
     'myeloid_epithelial_mixed_rate','median_doublet_score','median_contamination'),
     names(retention)),drop=FALSE]
writecsv(mixed_tbl,'04_Identity_MixedSignal_and_Retention.csv')

# Per-cell DecontX markers generated by previous workflow (not rerunning DecontX).
ambient_markers_path <- file.path(QC,'03_Stability','05_Selected_State_Ambient_Markers.csv')
marker_ambient <- read_checked(ambient_markers_path, required=FALSE)
if (!is.null(marker_ambient)) {
  stopifnot(all(c('cell','gene','raw_count','decont_count') %in% colnames(marker_ambient)))
  dat <- merge(marker_ambient,
       mapping[,c('cell','cluster','doublet_class','contamination'),drop=FALSE],
       by='cell',all.x=FALSE,all.y=FALSE)
  relevant <- unique(unlist(key_genes))
  relevant <- intersect(c(relevant,'EPCAM','KRT8','KRT18','KRT19','GAST','SST','CD3D'),unique(dat$gene))
  dat <- dat[dat$gene %in% relevant & dat$cluster %in% selected,,drop=FALSE]
  dc <- list()
  for (sc in scenario_names) {
    keep <- names(scenarios[[sc]])[scenarios[[sc]]]
    z <- dat[dat$cell %in% keep,,drop=FALSE]
    for (cl in selected) for (g in relevant) {
      zz <- z[z$cluster==cl & z$gene==g,,drop=FALSE]
      if (!nrow(zz)) next
      dc[[length(dc)+1L]] <- data.frame(
        scenario=sc,cluster=cl,state=states[[cl]],gene=g,n_cells=nrow(zz),
        frac_detected_raw=mean(zz$raw_count > 0),
        frac_detected_corrected=mean(zz$decont_count > 0),
        total_raw_UMI=sum(zz$raw_count),total_corrected_UMI=sum(zz$decont_count),
        fraction_UMI_removed=safe_ratio(sum(zz$raw_count-zz$decont_count),sum(zz$raw_count)))
    }
  }
  if(length(dc)) writecsv(do.call(rbind,dc),'05_DecontX_Marker_RawVsCorrected_ByState.csv')
  rm(dat,marker_ambient);gc(verbose=FALSE)
} else logmsg('Missing 05_Selected_State_Ambient_Markers.csv; corrected marker evidence not assessed')

# ---- 3: re-clustering concordance (independent ARI formula, no mclust) ------
logmsg('CHECK3: original 17-cluster vs QC reclustering concordance')
sc_files <- c(
  exclude_predicted_doublets=file.path(QC,'03_Stability','05_Recluster_exclude_predicted_doublets.rds'),
  exclude_doublets_ambient_le_0.2=file.path(QC,'03_Stability','05_Recluster_exclude_doublets_ambient_le_0.2.rds'))
concord <- list(); cluster_cons <- list()
for (sc in names(sc_files)) {
  if (!file.exists(sc_files[[sc]])) {logmsg('Missing recluster RDS: ',sc_files[[sc]]);next}
  x <- readRDS(sc_files[[sc]])
  a <- x$cell_cluster_assignments
  stopifnot(all(c('cell','original_cluster','new_cluster') %in% names(a)),!anyDuplicated(a$cell))
  a$original_cluster <- as.character(a$original_cluster)
  a$new_cluster <- as.character(a$new_cluster)
  ar <- ari_calc(a$original_cluster,a$new_cluster)
  tab <- table(a$original_cluster,a$new_cluster)
  cross <- as.data.frame(tab,stringsAsFactors=FALSE)
  names(cross) <- c('original_cluster','new_cluster','n')
  cross <- cross[cross$n>0,,drop=FALSE]
  cross$scenario <- sc
  writecsv(cross,paste0('06_Cluster_CrossTab_',sc,'.csv'))
  for (cl in selected) {
    row <- if (cl %in% rownames(tab)) tab[cl,,drop=TRUE] else numeric()
    cluster_cons[[length(cluster_cons)+1L]] <- data.frame(
      scenario=sc,cluster=cl,state=states[[cl]],n_retained=sum(row),
      n_new_clusters_contributing=sum(row>0),dominant_new_cluster=if(length(row)) names(which.max(row)) else NA_character_,
      fraction_in_dominant_new_cluster=if(sum(row)>0) max(row)/sum(row) else NA_real_)
  }
  concord[[length(concord)+1L]] <- data.frame(
    scenario=sc, n_cells=nrow(a), original_clusters=length(unique(a$original_cluster)),
    reclustered_clusters=length(unique(a$new_cluster)), ARI=ar,
    original_cells_matched=sum(a$cell %in% mapping$cell))
  rm(x);gc(verbose=FALSE)
}
if(length(concord)) {
  concord <- do.call(rbind,concord)
  writecsv(concord,'07_AdjustedRandIndex_ClusterStability.csv')
  writecsv(do.call(rbind,cluster_cons),'08_FiveState_ClusterFragmentation.csv')
  print(concord,row.names=FALSE)
}

# ---- 4a: Fig4F GO robustness with new marker contrasts ---------------------
# Re-evaluates selected cluster vs all other candidate cells in each QC subset.
# Uses original cluster identities; DOES NOT relabel cells based on reclustering.
# Original plot terms are read from Fig4_Final, without changing source results.
logmsg('CHECK4A: Fig4F GO-BP robustness')
fig4final <- 'E:/download/GSE150290_audit/Fig4_Final'
selected_go_path <- file.path(fig4final,'Fig4D_Selected_GO_BP_Terms.csv')
original_go <- read_checked(selected_go_path,required=FALSE)
if (!is.null(original_go)) {
  writecsv(original_go,'09_Original_Fig4F_GO_SelectedTerms_Copy.csv')
  logmsg('Found original Fig4 GO selected-term reference: ',nrow(original_go),' rows')
} else logmsg('Original GO selected-term table not found; GO can still be rerun')

if (RUN_DE_AND_GO && requireNamespace('Seurat',quietly=TRUE) &&
    requireNamespace('clusterProfiler',quietly=TRUE) &&
    requireNamespace('org.Hs.eg.db',quietly=TRUE)) {
  universe <- rownames(ct)[as.numeric(Matrix::rowMeans(ct > 0)) >= .01]
  universe <- intersect(universe,AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype='SYMBOL'))
  logmsg('GO background genes=',length(universe))
  go_list <- list(); msummary <- list()
  for (sc in setdiff(scenario_names,'baseline')) {
    gop <- file.path(OUT,paste0('10_GO_',sc,'.csv'))
    mkp <- file.path(OUT,paste0('10_DE_Markers_',sc,'.csv'))
    if(file.exists(gop) && file.exists(mkp)) {
      logmsg('CACHED GO results: ',sc)
      cached <- tryCatch(utils::read.csv(gop,stringsAsFactors=FALSE), error=function(e) NULL)
      if (!is.null(cached)) { go_list[[sc]] <- cached; next }
    }
    keep <- names(scenarios[[sc]])[scenarios[[sc]]]
    subset_counts <- ct[,keep,drop=FALSE]
    meta <- data.frame(original_cluster=mapping$cluster[match(keep,mapping$cell)],row.names=keep)
    s <- SeuratObject::CreateSeuratObject(counts=subset_counts,meta.data=meta,
                                          min.cells=0,min.features=0)
    s <- Seurat::NormalizeData(s,verbose=FALSE)
    s <- SeuratObject::`Idents<-`(s, value=as.character(s$original_cluster))
    markers <- list(); go <- list()
    for (cl in selected) {
      logmsg('DE/GO ',sc,' cluster ',cl)
      if(sum(SeuratObject::Idents(s) == cl) < 10L) {logmsg('SKIP low cell number: ',cl);next}
      res <- tryCatch(
        Seurat::FindMarkers(s,ident.1=cl,only.pos=TRUE, min.pct=.20,
                            logfc.threshold=.25,test.use='wilcox',verbose=FALSE),
        error=function(e) {logmsg('FindMarkers error: ',conditionMessage(e));NULL})
      if(is.null(res) || !nrow(res)) next
      res$gene <- rownames(res);res$cluster <- cl;res$state <- states[[cl]]
      res$scenario <- sc
      markers[[cl]] <- res
      fc_col <- if ('avg_log2FC' %in% names(res)) 'avg_log2FC' else 'avg_logFC'
      good <- !is.na(res$p_val_adj) & res$p_val_adj < .05 &
              res[[fc_col]] > .25 & res$pct.1 >= .25
      markers_set <- intersect(unique(res$gene[good]),universe)
      msummary[[length(msummary)+1L]] <- data.frame(
        scenario=sc,cluster=cl,state=states[[cl]],
        selected_marker_count=length(markers_set),n_state_cells=sum(SeuratObject::Idents(s) == cl))
      if(length(markers_set) < 10) {logmsg('GO SKIP cluster ',cl,': <10 input markers');next}
      eg <- tryCatch(clusterProfiler::enrichGO(gene=markers_set,
          OrgDb=org.Hs.eg.db::org.Hs.eg.db,keyType='SYMBOL',universe=universe,
          ont='BP',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,
          minGSSize=10,maxGSSize=500),
          error=function(e) {logmsg('GO error: ',conditionMessage(e));NULL})
      if(is.null(eg)) next
      gg <- as.data.frame(eg)
      if(nrow(gg)) {
        gg$cluster <- cl;gg$state <- states[[cl]];gg$scenario <- sc
        go[[cl]] <- gg
      }
    }
    if(length(markers)) writecsv(do.call(rbind,markers),paste0('10_DE_Markers_',sc,'.csv'))
    if(length(go)) {
      gout <- do.call(rbind,go)
      writecsv(gout,paste0('10_GO_',sc,'.csv'))
      go_list[[sc]] <- gout
    } else {writecsv(data.frame(),paste0('10_GO_',sc,'.csv'))}
    rm(s,subset_counts,markers,go);gc(verbose=FALSE)
  }
  if(length(msummary)) writecsv(do.call(rbind,msummary),'11_GO_Marker_Input_Summary.csv')
  if (!is.null(original_go) && all(c('state','ID') %in% names(original_go)) && length(go_list)) {
    # Track each originally displayed term, not just terms that remain significant.
    rows <- list()
    old <- unique(original_go[,intersect(c('state','ID','Description','p.adjust','Count'),names(original_go)),drop=FALSE])
    if(!('Description' %in% names(old))) old$Description <- NA_character_
    for(sc in names(go_list)) {
      new <- go_list[[sc]]
      if(!all(c('state','ID','p.adjust','Count') %in% names(new))) next
      q <- merge(old,new[,c('state','ID','p.adjust','Count')],by=c('state','ID'),all.x=TRUE,
                 suffixes=c('_original','_postQC'))
      q$scenario <- sc
      q$still_significant <- !is.na(q$p.adjust_postQC) & q$p.adjust_postQC < .05 &
                              q$Count_postQC >= 3
      rows[[length(rows)+1L]] <- q
    }
    if(length(rows)) writecsv(do.call(rbind,rows),'12_Fig4F_OriginalTerms_Significance_AfterQC.csv')
  }
} else {
  logmsg('GO not rerun: disabled or missing Seurat, clusterProfiler, org.Hs.eg.db')
}

# ---- 4b: Fig4D/E Slingshot robustness: coverage + optional independent refit
logmsg('CHECK4B: Fig4D/E Slingshot coverage and optional refit')
trajpath <- file.path(fig4final,'Fig4C_Trajectory2495_Slingshot.rds')
if (file.exists(trajpath) && requireNamespace('slingshot',quietly=TRUE) &&
    requireNamespace('SingleCellExperiment',quietly=TRUE)) {
  tr <- readRDS(trajpath)
  trcells <- colnames(tr)
  logmsg('Original trajectory object cells=',length(trcells),
         '; of 3101 candidates matched=',sum(trcells %in% mapping$cell))
  lc <- slingshot::slingLineages(tr)
  logmsg('Original saved Slingshot lineages=',length(lc))
  sp <- tryCatch(slingshot::slingParams(tr),error=function(e) NULL)
  capture.output(str(sp, max.level=3),file=file.path(OUT,'13_Slingshot_Original_Parameters.txt'))
  if(length(lc)) writecsv(data.frame(lineage=seq_along(lc),
               ordered_clusters=vapply(lc,function(x) paste(x,collapse=' -> '),character(1))),
               '14_Slingshot_Original_LineagePaths.csv')
  ptime <- tryCatch(slingshot::slingPseudotime(tr,na=TRUE),error=function(e) NULL)
  cover <- list()
  if(!is.null(ptime)) {
    for(sc in scenario_names) {
      mask <- as.logical(scenarios[[sc]][trcells])
      mask[is.na(mask)] <- FALSE
      for(j in seq_len(ncol(ptime))) {
        assigned <- is.finite(ptime[,j])
        cover[[length(cover)+1L]] <- data.frame(
          scenario=sc,lineage=j,n_original_assigned=sum(assigned),
          n_assigned_retained=sum(assigned & mask),
          assigned_retention=safe_ratio(sum(assigned & mask),sum(assigned)),
          median_pseudotime_retained=safe_med(ptime[assigned & mask,j]))
      }
    }
  }
  if(length(cover)) writecsv(do.call(rbind,cover),'15_Slingshot_LineageCoverage_PostQC.csv')
  rnames <- SingleCellExperiment::reducedDimNames(tr)
  logmsg('Trajectory reducedDims=',paste(rnames,collapse=','))
  roots <- if(length(lc)) unique(vapply(lc,function(z) as.character(z[1]),character(1))) else character()
  root_param <- if(!is.null(sp)) sp$start.clus else NULL
  root <- if(length(root_param)==1 && !is.na(root_param)) as.character(root_param) else
          if(length(roots)==1) roots else NULL
  root_origin <- if(length(root_param)==1 && !is.na(root_param)) 'saved_slingParams' else
                   if(length(roots)==1) 'inferred_from_first_elements_of_saved_lineages' else 'unresolved'
  writecsv(data.frame(root=if(length(root)) root else NA_character_,
         root_origin=root_origin,n_orig_cells=length(trcells),n_orig_lineages=length(lc)),
         '16_Slingshot_Root_Provenance.csv')
  if(RUN_SLINGSHOT_REFIT && length(root)==1 && length(rnames)) {
    # Refit using the original saved UMAP coordinates and cluster assignments.
    rd <- if(any(grepl('umap',rnames,ignore.case=TRUE))) {
      rnames[which(grepl('umap',rnames,ignore.case=TRUE))[1]]
    } else rnames[1]
    embedding <- SingleCellExperiment::reducedDim(tr,rd)
    labels <- tryCatch(slingshot::slingClusterLabels(tr),error=function(e) NULL)
    if(!is.null(labels)) {
      if(is.matrix(labels) && nrow(labels)==length(trcells)) {
        lab <- labels
      } else if(length(labels)==length(trcells)) {
        lab <- as.character(labels)
      } else lab <- NULL
      if(!is.null(lab)) for(sc in setdiff(scenario_names,'baseline')) {
        paths_file <- file.path(OUT,paste0('17_Slingshot_RefitPaths_',sc,'.csv'))
        if(file.exists(paths_file)) {logmsg('CACHED trajectory refit: ',sc);next}
        mask <- as.logical(scenarios[[sc]][trcells]);mask[is.na(mask)] <- FALSE
        if(sum(mask)<100) {logmsg('Slingshot refit SKIP too few retained cells: ',sc);next}
        x <- embedding[mask,,drop=FALSE]
        lsub <- if(is.matrix(lab)) lab[mask,,drop=FALSE] else lab[mask]
        logmsg('Slingshot refit ',sc,' kept=',nrow(x),' root=',root,
               ' (',root_origin,')')
        refit <- tryCatch(slingshot::slingshot(x,clusterLabels=lsub,
                    start.clus=root),
                error=function(e) {logmsg('Slingshot refit failed: ',conditionMessage(e));NULL})
        if(is.null(refit)) next
        newl <- slingshot::slingLineages(refit)
        table_out <- data.frame(scenario=sc,lineage=seq_along(newl),
          ordered_clusters=vapply(newl,function(z) paste(z,collapse=' -> '),character(1)),
          n_input_cells=sum(mask),root_provenance=root_origin,
          reducedDim=rd)
        writecsv(table_out,paste0('17_Slingshot_RefitPaths_',sc,'.csv'))
        best <- lapply(seq_along(lc),function(i) {
          oldseq <- as.character(lc[[i]])
          overlaps <- vapply(newl,function(newseq) {
             z <- as.character(newseq)
             safe_ratio(length(intersect(oldseq,z)),length(union(oldseq,z)))
          },numeric(1))
          data.frame(scenario=sc,old_lineage=i,best_refit_lineage=if(length(overlaps)) which.max(overlaps) else NA_integer_,
             best_cluster_path_Jaccard=if(length(overlaps)) max(overlaps) else NA_real_)
        })
        if(length(best)) writecsv(do.call(rbind,best),paste0('18_Slingshot_LineageMatching_',sc,'.csv'))
        logmsg('Slingshot refit success ',sc,' new lineages=',length(newl))
        rm(refit);gc(verbose=FALSE)
      }
    } else logmsg('Cannot reconstruct saved Slingshot cluster assignments: refit skipped')
  } else logmsg('Slingshot refit skipped: missing root/reduction or disabled; coverage still saved')
} else logmsg('Slingshot trajectory file/package not available; cannot assess trajectory')

# ---- Verdict evidence table: requires human review, no automatic PASS ------
checks <- data.frame(
  Fig4_panel=c('B/C','A/B/C','A','F','D/E'),
  analysis=c('Marker persistence and enrichment after doublet/ambient filtering',
             'APOE+/TREM2+ residual myeloid/epithelial and tissue sample support',
             '17-cluster partition consistency based on matched filtered cells',
             'GO Biological Process re-enrichment using filtered state markers',
             'Slingshot lineage coverage and reference-guided independent refitting'),
  expected_outputs=c('01_, 02_, 05_','03_, 04_','06_, 07_, 08_','10_, 11_, 12_','13_ through 18_'),
  interpretation_status='REVIEW_RESULTS_NOT_AUTOMATICALLY_VALIDATED')
writecsv(checks,'99_Fig4_Interpretation_Checklist.csv')
logmsg('FINISHED | Outputs=',OUT)
cat('\nFIG4 ROBUSTNESS FINISHED. Review CSVs and log in:\n',OUT,'\n')

}

# ---- legacy_07_Fig4_Slingshot_Checks : Slingshot and corrected markers verification ----
legacy_07_Fig4_Slingshot_Checks <- function() {
# ============================================================================
# PLOS ONE gastric scRNA-seq | DecontX numerical audit + Slingshot provenance
# Does NOT overwrite original Fig1-Fig4 objects or saved QC results.
# Run while getwd() == .../Fig1_Original81779_Revised:
#   source(file.choose(), encoding = "UTF-8")
# ============================================================================
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
stopifnot(basename(ROOT) == "Fig1_Original81779_Revised")
QC <- file.path(ROOT, "QC_Methods_Audit_20261021", "FourModules_v1")
OUT <- file.path(QC, "06_DecontX_Slingshot_Verification_v1")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TR_PATH <- "E:/download/GSE150290_audit/Fig4_Final/Fig4C_Trajectory2495_Slingshot.rds"
MAC_PATH <- file.path(ROOT, "Fig4_Macrophage_Reanalysis", "Phase6_HighResolution",
                      "07_HighRes_Annotated_Candidates.rds")
REFIT_TRAJECTORIES <- TRUE    # Set FALSE to verify only; original trajectories untouched.
SEED <- 20261021L
set.seed(SEED)

logfile <- file.path(OUT, "00_Verification_Log.txt")
logmsg <- function(...) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               paste0(..., collapse = ""), sep = " | ")
  cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE)
}
save_csv <- function(x, fn) {
  utils::write.csv(x, file = file.path(OUT, fn), row.names = FALSE, na = "")
}
read_csv <- function(f) {
  if (!file.exists(f)) stop("Required file missing: ", f)
  utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
fmt <- function(x) {
  if (is.null(x) || !length(x)) return("NULL")
  paste(capture.output(str(x, max.level = 1L)), collapse = " ")
}
scalar_string <- function(x) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(NA_character_)
  paste(as.character(x), collapse = ";")
}
safe_ratio <- function(x, y) if (is.finite(y) && y > 0) x / y else NA_real_
if (!requireNamespace("Matrix", quietly = TRUE) ||
    !requireNamespace("SeuratObject", quietly = TRUE) ||
    !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
    !requireNamespace("slingshot", quietly = TRUE)) {
  stop("Requires Matrix, SeuratObject, SingleCellExperiment, and slingshot")
}

mapping <- read_csv(file.path(QC, "03_Stability", "01_All3101_Audit_Mapping.csv"))
stopifnot(nrow(mapping) == 3101L, !anyDuplicated(mapping$cell),
          all(c("cell", "cluster", "doublet_class", "contamination") %in% names(mapping)),
          !anyNA(mapping$doublet_class), !anyNA(mapping$contamination))
mapping$cluster <- as.character(mapping$cluster)
states <- c("0" = "FCN1+ monocyte-like", "1" = "SPP1+ macrophages",
            "2" = "HLA-II+/GAST+ ambiguous", "3" = "C1Q+ macrophages",
            "7" = "APOE+/TREM2+ macrophages", "11" = "IL1B+ inflammatory myeloid")
selected_clusters <- names(states)
scenarios <- list(
  baseline = setNames(rep(TRUE, nrow(mapping)), mapping$cell),
  exclude_predicted_doublets = setNames(mapping$doublet_class != "doublet", mapping$cell),
  exclude_doublets_ambient_le_0.2 = setNames(
    mapping$doublet_class != "doublet" & mapping$contamination <= 0.2,
    mapping$cell)
)

# --------------------------------------------------------------------------
# A. DecontX: corrected counts are numerical posterior estimates, not integer
# counts. " > 0 " is inappropriate when corrected values are tiny residuals.
# We report 1e-6 (numerical), >=0.5 (primary operational cutoff), >=1 (strict)
# plus summed count removal. These are reporting cutoffs, not DecontX defaults.
# --------------------------------------------------------------------------
logmsg("START A: DecontX corrected-marker numerical verification")
if (!file.exists(MAC_PATH)) stop("Missing macrophage candidate object: ", MAC_PATH)
mac <- readRDS(MAC_PATH)
raw_counts <- SeuratObject::GetAssayData(mac, assay = "RNA", layer = "counts")
stopifnot(ncol(mac) == 3101L, all(mapping$cell %in% colnames(raw_counts)))
marker_genes <- c("APOE", "TREM2", "APOC1", "C1QA", "C1QB", "C1QC", "SPP1",
                  "LST1", "TYROBP", "FCER1G", "CD68", "CSF1R", "EPCAM",
                  "KRT8", "KRT18", "KRT19", "PIGR", "PHGR1", "GAST", "SST",
                  "CHGA", "CHGB", "CD3D", "CD3E")
marker_genes <- intersect(marker_genes, rownames(raw_counts))
ambient_dir <- file.path(QC, "02_Ambient")
ambfiles <- list.files(ambient_dir, pattern = "^ambient_.*\\.rds$", full.names = TRUE)
if (length(ambfiles) != 32L) logmsg("WARNING: expected 32 ambient RDS, found ", length(ambfiles))
if (!length(ambfiles)) stop("No saved ambient_[sample].rds files found")

per_cell <- vector("list", length(ambfiles))
manifest <- vector("list", length(ambfiles))
for (i in seq_along(ambfiles)) {
  f <- ambfiles[i]
  logmsg("DecontX [", i, "/", length(ambfiles), "] ", basename(f))
  x <- readRDS(f)
  sm <- unique(as.character(x$metrics$sample))
  if (length(sm) != 1) stop("Sample identifier ambiguous in ", f)
  cleaned <- x$corrected_counts
  if (is.null(cleaned)) stop("Missing corrected_counts in ", f,
                             ". Re-run original ambient module with SAVE_AMBIENT_CORRECTED_COUNTS=TRUE")
  scells <- intersect(colnames(cleaned), mapping$cell)
  sample_cells <- intersect(as.character(x$metrics$cell), mapping$cell)
  if (!setequal(scells, sample_cells)) stop("Ambient metrics and corrected matrix cell names differ: ", sm)
  genes <- intersect(marker_genes, rownames(cleaned))
  if (!setequal(genes, marker_genes)) {
    logmsg("WARNING: absent genes in ", sm, ": ", paste(setdiff(marker_genes, genes), collapse = ","))
  }
  manifest[[i]] <- data.frame(sample = sm, n_total = ncol(cleaned),
                               n_candidate = length(scells), n_genes = length(genes),
                               file = basename(f))
  if (!length(scells) || !length(genes)) {
    rm(x, cleaned); gc(verbose = FALSE); next
  }
  before <- as.matrix(raw_counts[genes, scells, drop = FALSE])
  after <- as.matrix(cleaned[genes, scells, drop = FALSE])
  if (any(!is.finite(after))) stop("Nonfinite corrected count in sample ", sm)
  if (any(after < -1e-8)) stop("Negative corrected count in sample ", sm)
  ix <- which(abs(before - round(before)) > 1e-8)
  if (length(ix)) logmsg("WARNING: raw reference counts noninteger in ", sm)
  one <- data.frame(
    cell = rep(scells, each = length(genes)),
    gene = rep(genes, times = length(scells)),
    raw_count = as.numeric(before),
    corrected_count = as.numeric(after),
    sample = sm,
    stringsAsFactors = FALSE
  )
  per_cell[[i]] <- one
  rm(x, cleaned, before, after, one); gc(verbose = FALSE)
}
save_csv(do.call(rbind, manifest), "01_DecontX_Sample_Coverage.csv")
per_cell <- Filter(Negate(is.null), per_cell)
if (!length(per_cell)) stop("No candidate cells found in corrected per-sample matrices")
dat <- do.call(rbind, per_cell)
rm(per_cell); gc(verbose = FALSE)
if (length(unique(dat$cell)) != 3101L) {
  stop("Incomplete corrected coverage: ", length(unique(dat$cell)), "/3101 cells")
}
if (anyDuplicated(dat[, c("cell", "gene")])) stop("Duplicate cell-gene corrected results")
dat$cluster <- mapping$cluster[match(dat$cell, mapping$cell)]
dat$state <- unname(states[dat$cluster])
dat$doublet_class <- mapping$doublet_class[match(dat$cell, mapping$cell)]
dat$contamination <- mapping$contamination[match(dat$cell, mapping$cell)]
# Audit existing candidate counts against DecontX posterior estimates.
raw_mismatch <- aggregate(abs(raw_count - round(raw_count)) ~ sample, dat,
                          function(v) sum(v > 1e-8))
save_csv(raw_mismatch, "02_Input_Count_Integer_Check.csv")
# Keep a small cell-by-gene audit table for all focus clusters.
focus_dat <- dat[dat$cluster %in% selected_clusters, , drop = FALSE]
save_csv(focus_dat, "03_Focus_Clusters_Corrected_Marker_PerCell.csv")

sumrows <- list()
for (sc in names(scenarios)) {
  cells <- names(scenarios[[sc]])[scenarios[[sc]]]
  for (cl in selected_clusters) {
    z <- dat[dat$cluster == cl & dat$cell %in% cells, , drop = FALSE]
    if (!nrow(z)) next
    for (gn in unique(z$gene)) {
      xx <- z[z$gene == gn, , drop = FALSE]
      n <- nrow(xx)
      sumrows[[length(sumrows) + 1L]] <- data.frame(
        scenario = sc, cluster = cl, state = states[[cl]], gene = gn,
        n_cells = n, total_raw = sum(xx$raw_count),
        total_corrected = sum(xx$corrected_count),
        fraction_count_removed = safe_ratio(
          sum(xx$raw_count) - sum(xx$corrected_count), sum(xx$raw_count)),
        raw_detection = mean(xx$raw_count >= 1),
        corrected_detection_gt0 = mean(xx$corrected_count > 0),
        corrected_detection_gt1e_6 = mean(xx$corrected_count > 1e-6),
        corrected_detection_ge0_5 = mean(xx$corrected_count >= 0.5),
        corrected_detection_ge1 = mean(xx$corrected_count >= 1),
        max_corrected_count = max(xx$corrected_count),
        p95_corrected_count = as.numeric(quantile(xx$corrected_count, 0.95)),
        stringsAsFactors = FALSE
      )
    }
  }
}
decont_summary <- do.call(rbind, sumrows)
save_csv(decont_summary, "04_DecontX_Marker_Threshold_Sensitivity.csv")

# Raw-vs-corrected cell signatures: >=2 marker genes, assess cluster7 and 2.
sets <- list(
  myeloid = c("LST1", "TYROBP", "FCER1G", "CD68", "CSF1R"),
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "PIGR", "PHGR1"),
  endocrine = c("GAST", "SST", "CHGA", "CHGB"))
combined <- list()
for (cl in c("2", "7")) {
  cells <- unique(dat$cell[dat$cluster == cl])
  for (sc in names(scenarios)) {
    cells <- intersect(unique(dat$cell[dat$cluster == cl]), names(scenarios[[sc]])[scenarios[[sc]]])
    if (!length(cells)) next
    d <- dat[dat$cell %in% cells, , drop = FALSE]
    u <- unique(d$cell)
    flags <- data.frame(cell = u, stringsAsFactors = FALSE)
    for (typ in names(sets)) {
      genes <- intersect(sets[[typ]], unique(d$gene))
      dd <- d[d$gene %in% genes, , drop = FALSE]
      rr <- table(factor(dd$cell[dd$raw_count >= 1], levels = u))
      cc <- table(factor(dd$cell[dd$corrected_count >= 0.5], levels = u))
      flags[[paste0(typ, "_raw")]] <- as.integer(rr) >= 2
      flags[[paste0(typ, "_corrected")]] <- as.integer(cc) >= 2
    }
    combined[[length(combined)+1L]] <- data.frame(
      cluster=cl, state=states[[cl]], scenario=sc, n_cells=length(u),
      n_myeloid_raw=sum(flags$myeloid_raw),
      n_myeloid_corrected=sum(flags$myeloid_corrected),
      n_mixed_myeloid_epithelial_raw=sum(flags$myeloid_raw & flags$epithelial_raw),
      n_mixed_myeloid_epithelial_corrected=sum(flags$myeloid_corrected & flags$epithelial_corrected),
      n_mixed_myeloid_endocrine_raw=sum(flags$myeloid_raw & flags$endocrine_raw),
      n_mixed_myeloid_endocrine_corrected=sum(flags$myeloid_corrected & flags$endocrine_corrected)
    )
  }
}
save_csv(do.call(rbind, combined), "05_Cluster2_7_MixedLineage_BeforeAfter.csv")
logmsg("A finished. Corrected-count positivity re-evaluated at >1e-6, >=0.5, >=1")
rm(dat, focus_dat, raw_counts, mac); gc(verbose = FALSE)

# --------------------------------------------------------------------------
# B. Slingshot: audit *the saved fitted object* rather than guessing its root.
# 2495 input cells; check inclusion/exclusion, stored start.given, fitted data
# coordinates, actual clusterLabels, 9 lineages, and optional matched refits.
# --------------------------------------------------------------------------
logmsg("START B: Slingshot provenance")
if (!file.exists(TR_PATH)) stop("Saved Slingshot file not found: ", TR_PATH)
tr <- readRDS(TR_PATH)
trcells <- colnames(tr)
stopifnot(length(trcells) == 2495L, !anyDuplicated(trcells), all(trcells %in% mapping$cell))
original_lineages <- slingshot::slingLineages(tr)
sp <- slingshot::slingParams(tr)
if (is.null(sp)) sp <- list()

params <- data.frame(
  parameter = names(sp),
  value = vapply(sp, function(v) paste(capture.output(str(v, max.level = 1)), collapse = " "), character(1)),
  stringsAsFactors = FALSE
)
save_csv(params, "10_Slingshot_Fitted_Params.csv")
writeLines(capture.output(str(sp, max.level = 3)), file.path(OUT, "10_Slingshot_Fitted_Params_Full.txt"))

start_given <- isTRUE(sp$start.given)
end_given <- if (!is.null(sp$end.given)) any(sp$end.given %in% TRUE) else FALSE
root_fit <- scalar_string(sp$start.clus)
firsts <- if (length(original_lineages))
  unique(vapply(original_lineages, function(z) as.character(z[1]), character(1))) else character()
root_prov <- if (start_given) "USER_SPECIFIED_IN_ORIGINAL_FIT" else
  if (identical(sp$start.given, FALSE)) "ALGORITHM_SELECTED_IN_ORIGINAL_FIT" else
    "NOT_VERIFIED_AS_USER_SPECIFIED"
if (is.na(root_fit) && length(firsts) == 1L) {
  root_fit <- firsts[1]
  root_prov <- "INFERRED_FROM_SAVED_PATH_ORDER__NOT_PROVEN_USER_SPECIFIED"
}
root_report <- data.frame(
  n_input_cells=length(trcells), n_total_candidates=nrow(mapping),
  n_saved_lineages=length(original_lineages), fitted_start=root_fit,
  start_given=if (is.null(sp$start.given)) NA else start_given,
  fitted_root_provenance=root_prov, fitted_end=scalar_string(sp$end.clus),
  end_explicit=if (is.null(sp$end.given)) NA else end_given,
  first_clusters=scalar_string(firsts),
  stringsAsFactors = FALSE
)
save_csv(root_report, "11_Slingshot_Root_and_Fit_Status.csv")

paths <- data.frame(
  lineage=seq_along(original_lineages),
  path=vapply(original_lineages, function(z) paste(z, collapse = " -> "), character(1)),
  first=vapply(original_lineages, function(z) as.character(z[1]), character(1)),
  terminal=vapply(original_lineages, function(z) as.character(z[length(z)]), character(1)),
  stringsAsFactors = FALSE
)
save_csv(paths, "12_Original_Nine_Lineage_Paths.csv")

mapping$in_trajectory <- mapping$cell %in% trcells
stages <- if ("stage_audit" %in% names(mapping)) "stage_audit" else NA_character_
samples <- if ("sample_audit" %in% names(mapping)) "sample_audit" else NA_character_
base <- data.frame(cluster=mapping$cluster, included=mapping$in_trajectory, stage=if(length(stages) && !is.na(stages)) mapping[[stages]] else "NA", sample=if(length(samples) && !is.na(samples)) mapping[[samples]] else "NA")
clusterinput <- as.data.frame(table(cluster=base$cluster, in_trajectory=base$included))
save_csv(clusterinput, "13_Trajectory_Input_Cluster_Coverage.csv")
stageinput <- as.data.frame(table(stage=base$stage, in_trajectory=base$included))
save_csv(stageinput, "14_Trajectory_Input_Stage_Coverage.csv")
smplinput <- as.data.frame(table(sample=base$sample, in_trajectory=base$included))
save_csv(smplinput, "15_Trajectory_Input_Sample_Coverage.csv")

# Identify original fitting space and ensure same coordinates are used for refit.
coords <- tryCatch(slingshot::slingReducedDim(tr), error = function(e) NULL)
if (is.null(coords)) stop("Could not obtain actual fitted coordinates via slingReducedDim(tr)")
coords <- as.matrix(coords)
if (nrow(coords) != length(trcells)) stop("Fitted embedding has incorrect cell count")
if (is.null(rownames(coords))) rownames(coords) <- trcells
if (!identical(rownames(coords), trcells)) {
  if (setequal(rownames(coords), trcells)) coords <- coords[trcells, , drop = FALSE]
  else stop("Fitted embedding row names inconsistent with trajectory cell barcodes")
}
rdnames <- SingleCellExperiment::reducedDimNames(tr)
coords_match <- lapply(rdnames, function(nm) {
  y <- as.matrix(SingleCellExperiment::reducedDim(tr, nm))
  same_shape <- identical(dim(coords), dim(y))
  data.frame(reduction=nm, n_dimensions=ncol(y),
             matches_saved_fitted_embedding=if (same_shape) isTRUE(all.equal(
               unname(coords), unname(y), tolerance=1e-7)) else FALSE,
             stringsAsFactors=FALSE)
})
save_csv(do.call(rbind, coords_match), "16_Fitted_Embedding_Identity.csv")

lab <- slingshot::slingClusterLabels(tr)
if (is.matrix(lab)) {
  if (nrow(lab) != length(trcells)) stop("Bad number of cluster-label rows")
  if (is.null(rownames(lab))) rownames(lab) <- trcells
  if (!identical(rownames(lab), trcells)) {
    if (setequal(rownames(lab), trcells)) lab <- lab[trcells, , drop=FALSE]
    else stop("Bad cluster-label rownames")
  }
  hard <- isTRUE(all.equal(as.numeric(rowSums(lab)),
                           rep(1, nrow(lab)), tolerance=1e-6)) &&
          all(abs(apply(lab,1,max)-1)<1e-6)
  hard_id <- if (hard) colnames(lab)[max.col(lab, ties.method="first")] else
    rep(NA_character_, length(trcells))
} else if (length(lab) == length(trcells)) {
  hard <- TRUE; hard_id <- as.character(lab)
  names(hard_id) <- trcells
} else {
  stop("Unrecognized original slingClusterLabels structure")
}
if (is.null(names(hard_id))) names(hard_id) <- trcells
if (!hard) logmsg("WARNING: Slingshot uses soft cluster assignment; refit preserves matrix weights")
cluster_comparison <- data.frame(
  cell=trcells, original_candidate_cluster=mapping$cluster[match(trcells, mapping$cell)],
  sling_cluster=hard_id, stringsAsFactors=FALSE)
save_csv(cluster_comparison, "17_Sling_versus_Recluster_Label_Mapping.csv")
if (hard) save_csv(as.data.frame(table(original_candidate_cluster=cluster_comparison$original_candidate_cluster,
                                        sling_cluster=cluster_comparison$sling_cluster)),
                   "18_CrossTab_Sling_Original_Clusters.csv")

# Refit using original fitted coordinates, the SAME saved cluster identities
# and reference-fitted starting orientation. No biologically justified root
# is asserted when start.given=FALSE or unavailable.
if (REFIT_TRAJECTORIES) {
  if (is.na(root_fit) || !nzchar(root_fit) || grepl(";", root_fit, fixed=TRUE)) {
    logmsg("REFIT NOT DONE: No single identifiable fitted root. See 11_* root status")
  } else {
    oldptime <- slingshot::slingPseudotime(tr, na=TRUE)
    stopifnot(nrow(oldptime)==length(trcells))
    for (sc in c("exclude_predicted_doublets", "exclude_doublets_ambient_le_0.2")) {
      cache <- file.path(OUT, paste0("20_Refit_", sc, ".rds"))
      if (file.exists(cache)) {
        logmsg("CACHED refit: ", sc)
        fr <- readRDS(cache)
      } else {
        mask <- as.logical(scenarios[[sc]][trcells]); mask[is.na(mask)] <- FALSE
        keep <- trcells[mask]
        lsub <- if (is.matrix(lab)) lab[keep, , drop=FALSE] else lab[mask]
        active_labels <- if (is.matrix(lsub)) colnames(lsub)[colSums(lsub)>0] else unique(as.character(lsub))
        if (!root_fit %in% active_labels) {
          logmsg("REFIT SKIPPED ", sc, ": reference root cluster absent after QC")
          next
        }
        args <- list(data=coords[keep, ,drop=FALSE], clusterLabels=lsub,
                     start.clus=root_fit)
        # Retain original specified terminal constraints only when explicitly given.
        if (end_given && !is.null(sp$end.clus) && all(sp$end.clus %in% active_labels))
          args$end.clus <- sp$end.clus
        carry <- c("dist.method","use.median","omega","omega_scale",
                   "shrink","extend","reweight","reassign","thresh","maxit",
                   "stretch","approx_points","smoother","shrink.method","allow.breaks")
        for (nm in carry) if (!is.null(sp[[nm]]) && length(sp[[nm]])==1 &&
                              !is.na(sp[[nm]]) && is.atomic(sp[[nm]])) args[[nm]] <- sp[[nm]]
        logmsg("Slingshot REFIT ", sc, " n=", length(keep), " root=", root_fit,
               " status=", root_prov, " fitted_dimensions=",ncol(coords))
        fit <- tryCatch(do.call(slingshot::slingshot, args),
                        error=function(e) {logmsg("FAILED refit: ", conditionMessage(e)); NULL})
        if (is.null(fit)) next
        fr <- list(scenario=sc, old_root=root_fit, root_provenance=root_prov,
                   n_input=length(keep),
                   pathways=slingshot::slingLineages(fit),
                   pseudo=slingshot::slingPseudotime(fit, na=TRUE),
                   cellnames=keep)
        saveRDS(fr, cache, compress="gzip")
        rm(fit); gc(verbose=FALSE)
      }
      new_paths <- fr$pathways
      summary_paths <- data.frame(
        scenario=sc, lineage=seq_along(new_paths),
        path=vapply(new_paths, function(z) paste(z, collapse=" -> "), character(1)),
        stringsAsFactors=FALSE)
      save_csv(summary_paths, paste0("21_Refit_LineagePaths_",sc,".csv"))
      matched <- list()
      for (i in seq_along(original_lineages)) {
        oldseq <- as.character(original_lineages[[i]])
        jaccard <- vapply(new_paths, function(z) {
          seq2 <- as.character(z)
          safe_ratio(length(intersect(oldseq,seq2)), length(union(oldseq,seq2)))
        }, numeric(1))
        if (!length(jaccard)) next
        j <- which.max(jaccard)
        common <- intersect(trcells, fr$cellnames)
        original_values <- oldptime[match(common,trcells),i]
        new_values <- fr$pseudo[match(common,fr$cellnames),j]
        valid <- is.finite(original_values) & is.finite(new_values)
        rho <- if (sum(valid) >= 30L && length(unique(original_values[valid])) >= 3L &&
                   length(unique(new_values[valid])) >= 3L) {
          suppressWarnings(cor(original_values[valid],new_values[valid],method="spearman"))
        } else NA_real_
        matched[[length(matched)+1L]] <- data.frame(
          scenario=sc, original_lineage=i, closest_refit_lineage=j,
          original_path=paste(oldseq,collapse=" -> "),
          refit_path=paste(new_paths[[j]],collapse=" -> "),
          identical_ordered_path=identical(oldseq,as.character(new_paths[[j]])),
          cluster_set_Jaccard=max(jaccard), n_shared_pseudotime_cells=sum(valid),
          spearman_pseudotime=rho, root_provenance=root_prov,
          stringsAsFactors=FALSE)
      }
      if(length(matched)) save_csv(do.call(rbind,matched),
                                  paste0("22_Refit_PathAndPseudotime_",sc,".csv"))
    }
  }
}
logmsg("FINISHED | All original sources preserved | output=", OUT)
cat("\nVERIFICATION COMPLETE. See:\n", OUT, "\n")

}

# ---- legacy_08_Fig5_AUCell_Plotting_v2 : Historical Fig5 A-F v2, differs from FINAL D/E ----
legacy_08_Fig5_AUCell_Plotting_v2 <- function() {
# ============================================================
# FIGURE 5 | Myeloid TME remodeling and candidate TF regulation
#
# Run in the existing Windows R session:
#   source(file.choose(), encoding = "UTF-8")
#
# Inputs: cached corrected Seurat object, cached 326x1529 AUCell
#         matrix, cached SCENIC regulons, existing `fit_D`.
# No SCENIC or AUCell recalculation; no single-cell p-values.
# A/B/C: TF states and stages; D/E: tissue-stage myeloid remodeling;
# F: exploratory network of 50 candidate TF regulons.
# ============================================================

local({
  needed <- c("ggplot2", "SeuratObject", "slingshot", "Matrix", "igraph")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Please install missing packages: ", paste(missing, collapse = ", "))

  RUN_DIR <- "E:/download/GSE150290_audit/Fig4_Final/Fig4_Manuscript_Rebuild_v1"
  TF_DIR <- file.path(RUN_DIR, "08_TF_SixStates")
  OUT_DIR <- file.path(TF_DIR, "07_Fig5_ABCDEF_TME_v2_1529")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("\nFIG5 | TME remodeling A-F v2, audited 1529 cells\n")

  auc_path <- file.path(TF_DIR, "04_SCENIC_Regulon_AUCell_1529.rds")
  seu_path <- file.path(TF_DIR, "00_TF_Input_1529_DecontX_WithSixStates.rds")
  sets_path <- file.path(TF_DIR, "01_Reused_SCENIC_Regulons.rds")
  stopifnot(file.exists(auc_path), file.exists(seu_path), file.exists(sets_path))
  auc <- readRDS(auc_path)
  if (inherits(auc, "aucellResults")) {
    if (!requireNamespace("AUCell", quietly = TRUE)) stop("AUCell not installed")
    auc <- AUCell::getAUC(auc)
  }
  auc <- as.matrix(auc)
  seu <- readRDS(seu_path)
  stopifnot(inherits(seu, "Seurat"), ncol(seu) == 1529L,
            ncol(auc) == 1529L, setequal(colnames(auc), colnames(seu)),
            !anyNA(auc), all(is.finite(auc)))
  cells <- colnames(seu)
  auc <- auc[, cells, drop = FALSE]

  meta <- seu@meta.data[cells, , drop = FALSE]
  stopifnot(all(c("fig4_state_display", "stage", "sample_uid") %in% colnames(meta)))
  meta$cell <- cells
  meta$state <- as.character(meta$fig4_state_display)
  meta$stage <- as.character(meta$stage)
  meta$sample_uid <- as.character(meta$sample_uid)
  stages <- c("NAG", "CAG", "IM", "GC")
  if (anyNA(meta$state) || anyNA(meta$stage) || anyNA(meta$sample_uid) ||
      !all(meta$stage %in% stages)) stop("Unexpected missing/invalid state, stage or sample metadata")

  states5 <- c("FCN1+ monocyte-like", "IL1B+ inflammatory myeloid",
               "SPP1+ macrophages", "APOE+/TREM2+ macrophages", "C1Q+ macrophages")
  stopifnot(all(states5 %in% unique(meta$state)))
  mixed <- setdiff(unique(meta$state), states5)
  if (length(mixed) != 1L || !grepl("HLA-II|GAST", mixed)) {
    stop("Expected one HLA-II/GAST mixed-marker state; found: ", paste(mixed, collapse = ", "))
  }
  states6 <- c(states5, mixed)
  stopifnot(length(states6) == 6L, all(table(factor(meta$state, levels = states6)) > 0L))
  short_state <- setNames(
    c("FCN1+", "IL1B+", "SPP1+", "APOE+/TREM2+", "C1Q+", "HLA-II+/GAST+ (mixed)"), states6)
  palette6 <- setNames(
    c("#B88238", "#9F77C9", "#C55286", "#4EAB8D", "#4795BA", "#A19591"), states6)
  palette_stage <- c("NAG" = "#4795BA", "CAG" = "#4EAB8D",
                     "IM" = "#B88238", "GC" = "#C55286")

  # GC technical samples from the same numerical donor are collapsed.
  # Other cohorts remain at the independently labeled sample level, unless
  # an authoritative donor mapping becomes available.
  sample_tag <- sub("^.*::", "", meta$sample_uid)
  if (any(meta$stage == "GC" & !grepl("^[0-9]+_t[0-9]+$", sample_tag))) {
    stop("GC donor identity cannot be parsed from sample_uid; do not guess donor identities")
  }
  meta$unit <- ifelse(meta$stage == "GC",
                      paste0("GC::", sub("_t[0-9]+$", "", sample_tag)),
                      paste0(meta$stage, "::", meta$sample_uid))
  meta$key <- paste(meta$unit, meta$state, sep = "||")
  unit_meta <- unique(meta[, c("unit", "stage"), drop = FALSE])
  key_meta <- unique(meta[, c("key", "unit", "stage", "state"), drop = FALSE])
  if (anyDuplicated(unit_meta$unit) || anyDuplicated(key_meta$key)) {
    stop("Nonunique biological unit identifier detected")
  }
  stage_audit <- as.data.frame(table(stage = factor(meta$stage, levels = stages),
                                     state = factor(meta$state, levels = states6)))
  utils::write.csv(stage_audit, file.path(OUT_DIR, "00_Stage_By_State_CellCount_Audit.csv"),
                   row.names = FALSE)
  unit_audit <- as.data.frame(table(stage = factor(unit_meta$stage, levels = stages)))
  utils::write.csv(unit_audit, file.path(OUT_DIR, "00_Independent_Unit_Counts_By_Stage.csv"),
                   row.names = FALSE)

  theme_sci <- function() {
    ggplot2::theme_classic(base_size = 11, base_family = "sans") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#58656D"),
        plot.caption = ggplot2::element_text(size = 8.5, colour = "#62717A", hjust = 0),
        axis.title = ggplot2::element_text(size = 10.5),
        axis.text = ggplot2::element_text(size = 9.5, colour = "#343D44"),
        strip.background = ggplot2::element_rect(fill = "#F0F3F4", colour = NA),
        strip.text = ggplot2::element_text(face = "bold", size = 10),
        legend.title = ggplot2::element_text(face = "bold", size = 9.5),
        plot.margin = ggplot2::margin(12, 15, 12, 12)
      )
  }
  save_fig <- function(p, label, w, h) {
    png <- file.path(OUT_DIR, paste0(label, ".png"))
    pdf <- file.path(OUT_DIR, paste0(label, ".pdf"))
    ggplot2::ggsave(png, p, width = w, height = h, units = "in", dpi = 500, bg = "white")
    ggplot2::ggsave(pdf, p, width = w, height = h, units = "in", bg = "white")
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggplot2::ggsave(file.path(OUT_DIR, paste0(label, ".svg")), p,
                      width = w, height = h, units = "in", device = "svg", bg = "white")
    } else {
      message("SVG skipped for ", label, ": install svglite to enable SVG")
    }
    cat("Saved ", label, "\n", sep = "")
  }
  zscore <- function(v) {
    s <- stats::sd(v)
    if (!is.finite(s) || s < 1e-12) return(rep(0, length(v)))
    (v - mean(v)) / s
  }
  group_means <- function(mat, groups) {
    stopifnot(nrow(mat) == length(groups), !anyNA(groups))
    sizes <- table(groups)
    sums <- rowsum(mat, group = groups, reorder = TRUE)
    out <- sweep(sums, 1L, as.numeric(sizes[rownames(sums)]), "/")
    list(means = out, n = as.integer(sizes[rownames(sums)]))
  }

  # One standard (non-extended if possible) regulon per prespecified TF.
  info <- data.frame(regulon = rownames(auc), stringsAsFactors = FALSE)
  info$TF <- sub("(_extended)? \\([0-9]+g\\)$", "", info$regulon)
  info$n_targets <- suppressWarnings(as.integer(
    sub("^.*\\(([0-9]+)g\\)$", "\\1", info$regulon)))
  info$extended <- grepl("_extended \\(", info$regulon)
  info <- info[is.finite(info$n_targets) & info$n_targets >= 20L &
                 info$n_targets <= 1000L, , drop = FALSE]
  tf_pref <- c("MAF", "JUNB", "STAT1", "REL", "CEBPB", "FOSL2",
               "IRF1", "TFEB", "XBP1", "NFE2L2", "CREM", "SPI1")
  selected <- do.call(rbind, lapply(tf_pref, function(nm) {
    z <- info[info$TF == nm, , drop = FALSE]
    if (!nrow(z)) return(NULL)
    z[order(z$extended, z$n_targets), , drop = FALSE][1L, , drop = FALSE]
  }))
  if (is.null(selected) || nrow(selected) < 8L) {
    stop("Too few prespecified TF regulons with 20-1000 targets")
  }
  tf4 <- c("REL", "MAF", "STAT1", "CEBPB")
  if (!all(c(tf4, "TFEB") %in% selected$TF)) {
    stop("Required REL / MAF / STAT1 / CEBPB / TFEB regulon missing")
  }
  utils::write.csv(selected, file.path(OUT_DIR, "00_Selected_Regulons.csv"), row.names = FALSE)

  # Cell-by-TF matrices, both AUCell and comparable within-regulon Z values.
  X <- t(auc[selected$regulon, cells, drop = FALSE])
  rownames(X) <- cells
  colnames(X) <- selected$TF
  Xz <- apply(X, 2L, zscore)
  Xz <- as.matrix(Xz)
  dimnames(Xz) <- dimnames(X)
  state_raw_mean <- group_means(X, meta$state)$means
  state_z_mean <- group_means(Xz, meta$state)$means
  Xr <- Xz - state_z_mean[meta$state, , drop = FALSE]

  # ==========================================================
  # A | Six-state TF-regulon activity heatmap (keep approved layout)
  # ==========================================================
  A <- as.data.frame(as.table(state_raw_mean), stringsAsFactors = FALSE)
  names(A) <- c("state", "TF", "mean_AUCell")
  A$state <- as.character(A$state)
  A$TF <- as.character(A$TF)
  A$z <- ave(A$mean_AUCell, A$TF, FUN = function(x) pmax(-2, pmin(2, zscore(x))))
  A$state <- factor(A$state, levels = states6)
  A$TF <- factor(A$TF, levels = rev(selected$TF))
  pA <- ggplot2::ggplot(A, ggplot2::aes(x = state, y = TF, fill = z)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.65, width = 0.98, height = 0.97) +
    ggplot2::scale_x_discrete(labels = short_state, drop = FALSE) +
    ggplot2::scale_fill_gradient2(low = "#4795BA", mid = "#F7F8F8",
                                  high = "#C55286", midpoint = 0,
                                  limits = c(-2, 2), name = "Within-regulon Z") +
    ggplot2::labs(title = "A  |  State-specific TF regulon landscape",
                  subtitle = "Six audited states across NAG, CAG, IM and GC",
                  x = NULL, y = NULL,
                  caption = "Descriptive cell means; HLA-II+/GAST+ is a mixed-marker state.") +
    theme_sci() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9))
  utils::write.csv(A, file.path(OUT_DIR, "Fig5A_SixState_Regulons.csv"), row.names = FALSE)
  save_fig(pA, "Fig5A_SixState_Regulon_Heatmap", 12, 6.5)

  # ==========================================================
  # B | Same approved UMAP TF panels, same Fig4D coordinates
  # ==========================================================
  if (!exists("fit_D", envir = .GlobalEnv, inherits = FALSE)) {
    stop("fit_D is not in your R session; please load 1529-cell Fig4D fit_D first")
  }
  fit <- get("fit_D", envir = .GlobalEnv)
  xy <- as.matrix(slingshot::slingReducedDim(fit))
  if (is.null(rownames(xy)) || ncol(xy) < 2L ||
      !setequal(rownames(xy), cells)) stop("fit_D coordinates do not match 1529 cells")
  xy <- xy[cells, 1:2, drop = FALSE]
  B <- do.call(rbind, lapply(tf4, function(tf) {
    data.frame(cell = cells, x = xy[, 1L], y = xy[, 2L], TF = tf,
               AUC = X[, tf], z = pmax(-2, pmin(2, Xz[, tf])),
               stringsAsFactors = FALSE)
  }))
  B$TF <- factor(B$TF, levels = tf4)
  B <- B[order(B$TF, B$z), , drop = FALSE]
  pB <- ggplot2::ggplot(B, ggplot2::aes(x = x, y = y, colour = z)) +
    ggplot2::geom_point(size = 0.65, alpha = 0.9, stroke = 0) +
    ggplot2::facet_wrap(~ TF, ncol = 2) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_gradient2(low = "#4795BA", mid = "#E9ECEB",
                                    high = "#C55286", midpoint = 0,
                                    limits = c(-2, 2), name = "Within-regulon Z") +
    ggplot2::labs(title = "B  |  Spatial distribution of TF activity",
                  subtitle = "Corrected 1529 cells; frozen Fig4D coordinates",
                  x = "UMAP 1", y = "UMAP 2",
                  caption = "Colours indicate inferred regulon activity, not TF mRNA expression.") +
    theme_sci() + ggplot2::theme(panel.spacing = grid::unit(1, "lines"))
  utils::write.csv(B, file.path(OUT_DIR, "Fig5B_TF_UMAP_Data.csv"), row.names = FALSE)
  save_fig(pB, "Fig5B_TF_Regulon_UMAP", 9.5, 8)

  # ==========================================================
  # C | NAG > CAG > IM > GC: stage changes at independent-unit level
  #
  # First center within the six states (Xr); then average by independent
  # biological unit. This reduces, but does not eliminate, differences
  # caused simply by changing six-state composition among stages.
  # ==========================================================
  unit_result <- group_means(Xr, meta$unit)
  unit_wide <- unit_result$means
  C <- as.data.frame(as.table(unit_wide), stringsAsFactors = FALSE)
  names(C) <- c("unit", "TF", "state_adjusted_Z")
  C$unit <- as.character(C$unit)
  C$TF <- as.character(C$TF)
  C$stage <- unit_meta$stage[match(C$unit, unit_meta$unit)]
  C$n_cells <- unit_result$n[match(C$unit, rownames(unit_wide))]
  C <- C[C$n_cells >= 5L, , drop = FALSE]
  C$stage <- factor(C$stage, levels = stages)
  C$TF <- factor(C$TF, levels = tf4)
  C4 <- C[!is.na(C$TF), , drop = FALSE]
  if (!nrow(C4) || length(unique(C4$stage)) < 3L) {
    stop("Insufficient stage support for Fig5C after minimum 5-cell biological-unit gate")
  }
  pC <- ggplot2::ggplot(C4, ggplot2::aes(x = stage, y = state_adjusted_Z)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#A7AFB5", linetype = "dashed") +
    ggplot2::geom_boxplot(width = 0.51, outlier.shape = NA,
                          colour = "#87949B", fill = "#E6EBED", alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(colour = stage), size = 2.0,
                        alpha = 0.92, position = ggplot2::position_jitter(width = 0.08, height = 0)) +
    ggplot2::facet_wrap(~ TF, ncol = 2, scales = "fixed") +
    ggplot2::scale_colour_manual(values = palette_stage, drop = FALSE) +
    ggplot2::labs(title = "C  |  TF regulon activity across disease stages",
                  subtitle = "Independent GC donors / other labeled samples; state-centered activity",
                  x = "Histological stage", y = "Mean state-adjusted AUCell Z",
                  caption = "Each point is one independent unit; distributions are descriptive (no single-cell p-values).") +
    theme_sci() + ggplot2::theme(legend.position = "none", panel.spacing = grid::unit(1, "lines"))
  utils::write.csv(C, file.path(OUT_DIR, "Fig5C_Stage_IndependentUnit_Activity.csv"), row.names = FALSE)
  save_fig(pC, "Fig5C_Stage_Regulon_DonorDots", 9.5, 8)

  # ==========================================================
  # D | Six-state composition across NAG, CAG, IM and GC
  # Within the selected myeloid-associated compartment, not whole-tissue
  # macrophage infiltration. Each bar is one donor / labeled sample.
  # ==========================================================
  state_count <- as.data.frame(table(unit = meta$unit, state = factor(meta$state, levels = states6)),
                               stringsAsFactors = FALSE)
  names(state_count) <- c("unit", "state", "n_cells")
  state_count$unit <- as.character(state_count$unit)
  state_count$state <- as.character(state_count$state)
  state_count$stage <- unit_meta$stage[match(state_count$unit, unit_meta$unit)]
  state_count$total_myeloid <- as.numeric(tapply(state_count$n_cells, state_count$unit, sum)[state_count$unit])
  state_count$fraction <- state_count$n_cells / state_count$total_myeloid
  D <- state_count[state_count$total_myeloid >= 5L, , drop = FALSE]
  if (!nrow(D) || length(unique(D$stage)) < 3L) {
    stop("Insufficient biological-unit support for Fig5D six-state composition")
  }
  unit_order <- unique(D[, c("unit", "stage", "total_myeloid")])
  unit_order <- unit_order[order(match(unit_order$stage, stages),
                                 -unit_order$total_myeloid, unit_order$unit), , drop = FALSE]
  D$unit <- factor(D$unit, levels = unit_order$unit)
  D$stage <- factor(D$stage, levels = stages)
  D$state <- factor(D$state, levels = states6)
  # Short labels preserve distinct sample/donor units without inventing IDs.
  shown_id <- setNames(paste0("U", seq_len(nrow(unit_order))), unit_order$unit)
  utils::write.csv(transform(unit_order, plot_label = shown_id[unit]),
                   file.path(OUT_DIR, "Fig5D_Unit_Label_Map.csv"), row.names = FALSE)
  pD <- ggplot2::ggplot(D, ggplot2::aes(x = unit, y = fraction, fill = state)) +
    ggplot2::geom_col(width = 0.88, colour = "white", linewidth = 0.13) +
    ggplot2::facet_grid(~ stage, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_discrete(labels = shown_id, drop = TRUE) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                                limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.025))) +
    ggplot2::scale_fill_manual(values = palette6, labels = short_state, name = "Six states", drop = FALSE) +
    ggplot2::labs(title = "D  |  Myeloid-state composition across tissue stages",
                  subtitle = "Each column denotes one independent donor / labeled sample",
                  x = NULL, y = "Within-compartment proportion",
                  caption = "Relative fractions among selected six states; not whole-tissue macrophage abundance. Units with <5 cells omitted.") +
    theme_sci() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 7.3),
                                 panel.spacing.x = grid::unit(0.9, "lines"),
                                 legend.position = "bottom",
                                 legend.text = ggplot2::element_text(size = 8)) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
  utils::write.csv(D, file.path(OUT_DIR, "Fig5D_SixState_Composition_ByUnit.csv"), row.names = FALSE)
  save_fig(pD, "Fig5D_Myeloid_State_Composition", 13.1, 6.1)

  # ==========================================================
  # E | Functional immune program changes across tissue stages
  # Each point is a donor / labeled sample (not an individual cell).
  # Unlike C, these are NOT residualized for cell-state composition;
  # they capture both changes in state mixture and within-state programs.
  # ==========================================================
  panels <- list(
    Inflammatory = c("IL1B", "TNF", "CCL3", "CCL4", "CXCL8", "NFKBIA", "PTGS2"),
    Antigen_presentation = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA"),
    Lysosomal_processing = c("CTSD", "CTSB", "CTSL", "LAMP1", "LAMP2", "LIPA"),
    Lipid_handling = c("APOE", "TREM2", "LPL", "GPNMB", "ABCA1", "APOC1")
  )
  counts <- tryCatch(
    SeuratObject::GetAssayData(seu, assay = "RNA", layer = "counts"),
    error = function(e) SeuratObject::GetAssayData(seu, assay = "RNA", slot = "counts")
  )
  stopifnot(setequal(colnames(counts), cells))
  counts <- counts[, cells, drop = FALSE]
  libsize <- Matrix::colSums(counts)
  stopifnot(all(is.finite(libsize)), all(libsize > 0))
  gene_audit <- do.call(rbind, lapply(names(panels), function(nm) {
    data.frame(program = nm, gene = panels[[nm]],
               detected = panels[[nm]] %in% rownames(counts))
  }))
  utils::write.csv(gene_audit, file.path(OUT_DIR, "Fig5E_Functional_Marker_Gene_Coverage.csv"),
                   row.names = FALSE)
  if (any(vapply(panels, function(z) sum(z %in% rownames(counts)) < 3L, logical(1)))) {
    stop("One or more functional gene sets have fewer than 3 genes in corrected RNA counts")
  }
  valid_genes <- unique(unlist(lapply(panels, intersect, y = rownames(counts))))
  gene_counts <- as.matrix(counts[valid_genes, , drop = FALSE])
  gene_logcpm <- log1p(t(t(gene_counts) / as.numeric(libsize)) * 1e4)
  module_score <- do.call(cbind, lapply(panels, function(z) {
    colMeans(gene_logcpm[intersect(z, rownames(gene_logcpm)), , drop = FALSE])
  }))
  rownames(module_score) <- cells
  colnames(module_score) <- names(panels)
  module_z <- apply(module_score, 2L, zscore)
  module_z <- as.matrix(module_z)
  dimnames(module_z) <- dimnames(module_score)
  module_units <- group_means(module_z, meta$unit)
  E <- as.data.frame(as.table(module_units$means), stringsAsFactors = FALSE)
  names(E) <- c("unit", "program", "mean_expression_Z")
  E$unit <- as.character(E$unit)
  E$program <- as.character(E$program)
  E$stage <- unit_meta$stage[match(E$unit, unit_meta$unit)]
  E$n_cells <- module_units$n[match(E$unit, rownames(module_units$means))]
  E <- E[E$n_cells >= 5L & is.finite(E$mean_expression_Z), , drop = FALSE]
  if (!nrow(E) || length(unique(E$stage)) < 3L) {
    stop("Insufficient biological-unit support for Fig5E functional stages")
  }
  program_label <- c(Inflammatory = "Inflammatory signaling",
                     Antigen_presentation = "Antigen presentation",
                     Lysosomal_processing = "Lysosomal processing",
                     Lipid_handling = "Lipid handling")
  E$program <- factor(E$program, levels = names(program_label), labels = unname(program_label))
  E$stage <- factor(E$stage, levels = stages)
  E_groups <- split(E, interaction(E$program, E$stage, drop = TRUE))
  E_summary <- do.call(rbind, lapply(E_groups, function(d) {
    v <- d$mean_expression_Z
    data.frame(program = as.character(d$program[1L]),
               stage = as.character(d$stage[1L]), n_units = length(v),
               median = stats::median(v),
               q25 = as.numeric(stats::quantile(v, 0.25)),
               q75 = as.numeric(stats::quantile(v, 0.75)))
  }))
  E_summary$program <- factor(E_summary$program, levels = levels(E$program))
  E_summary$stage <- factor(E_summary$stage, levels = stages)
  pE <- ggplot2::ggplot(E, ggplot2::aes(x = mean_expression_Z, y = stage)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BCC6CA", linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(colour = stage), alpha = 0.68, size = 2.1,
                        position = ggplot2::position_jitter(width = 0, height = 0.12)) +
    ggplot2::geom_segment(data = E_summary,
                          ggplot2::aes(x = q25, xend = q75, y = stage, yend = stage),
                          inherit.aes = FALSE, linewidth = 2.7, colour = "#485A63") +
    ggplot2::geom_point(data = E_summary, ggplot2::aes(x = median, y = stage),
                        inherit.aes = FALSE, shape = 21, size = 3.6,
                        stroke = 0.85, fill = "white", colour = "#31454F") +
    ggplot2::facet_wrap(~ program, ncol = 2, scales = "free_x") +
    ggplot2::scale_colour_manual(values = palette_stage, drop = FALSE) +
    ggplot2::labs(title = "E  |  Functional remodeling of the myeloid compartment",
                  subtitle = "Immune gene programs across NAG, CAG, IM and GC",
                  x = "Mean immune-module expression (cell-standardized Z)",
                  y = "Histological stage",
                  caption = "Faint points = independent donors / labeled samples; white circles = medians; bars = interquartile ranges. Descriptive only.") +
    theme_sci() + ggplot2::theme(panel.spacing = grid::unit(1.1, "lines"),
                                 legend.position = "none")
  utils::write.csv(E, file.path(OUT_DIR, "Fig5E_Functional_Modules_ByUnit.csv"), row.names = FALSE)
  utils::write.csv(E_summary, file.path(OUT_DIR, "Fig5E_Functional_Modules_ByStage_Summary.csv"), row.names = FALSE)
  save_fig(pE, "Fig5E_Myeloid_Immune_Programs_ByStage", 10.4, 8.0)

  # ==========================================================
  # F | Top 50 candidate TFs, ranked by state and stage variation.
  # Nodes: one preferred (non-extended if present) regulon per TF.
  # Edges: shared SCENIC target genes (NOT TF-TF regulation).
  # Node colour: exploratory stage of highest within-state-centered activity;
  # low-contrast or poorly supported stages receive neutral grey.
  # ==========================================================
  regulon_list <- readRDS(sets_path)
  if (!is.list(regulon_list) || is.null(names(regulon_list))) {
    stop("Saved SCENIC regulon network is not a named list")
  }
  as_genes <- function(x) {
    if (is.character(x)) return(unique(x[nzchar(x)]))
    if (is.numeric(x) && !is.null(names(x))) return(unique(names(x)))
    if (is.list(x) && !is.null(x$genes)) return(unique(as.character(x$genes)))
    if (is.list(x)) {
      z <- unlist(x, use.names = FALSE)
      if (is.character(z)) return(unique(z[nzchar(z)]))
    }
    character(0)
  }
  all_info <- info[info$regulon %in% names(regulon_list), , drop = FALSE]
  all_info <- all_info[order(all_info$TF, all_info$extended,
                             all_info$n_targets, all_info$regulon), , drop = FALSE]
  all_info <- all_info[!duplicated(all_info$TF), , drop = FALSE]
  all_targets <- lapply(all_info$regulon, function(r) as_genes(regulon_list[[r]]))
  good_targets <- lengths(all_targets) >= 20L
  all_info <- all_info[good_targets, , drop = FALSE]
  all_targets <- all_targets[good_targets]
  if (nrow(all_info) < 50L) {
    stop("F requires 50 distinct TFs with usable 20-1000-target regulons; available: ",
         nrow(all_info), ". Do not fabricate missing TFs.")
  }
  # All regulons were normalized independently across cells.
  candidate_auc <- auc[all_info$regulon, cells, drop = FALSE]
  candidate_z <- t(apply(candidate_auc, 1L, zscore))
  dimnames(candidate_z) <- dimnames(candidate_auc)
  candidate_cells <- t(candidate_z)
  colnames(candidate_cells) <- all_info$TF
  cell_state_means <- group_means(candidate_cells, meta$state)$means
  candidate_centered <- candidate_cells - cell_state_means[meta$state, , drop = FALSE]
  candidate_units <- group_means(candidate_centered, meta$unit)
  unit_ok <- candidate_units$n >= 5L
  candidate_unit_stages <- unit_meta$stage[match(rownames(candidate_units$means), unit_meta$unit)]
  stage_support <- vapply(stages, function(st) sum(unit_ok & candidate_unit_stages == st), integer(1))
  if (any(stage_support < 2L)) {
    stop("Insufficient stage representation for ranking 50 TFs: ",
         paste(paste(stages, stage_support, sep = "="), collapse = ", "))
  }
  unit_stage_matrix <- vapply(stages, function(st) {
    colMeans(candidate_units$means[unit_ok & candidate_unit_stages == st, , drop = FALSE])
  }, numeric(nrow(all_info)))
  rownames(unit_stage_matrix) <- all_info$TF
  colnames(unit_stage_matrix) <- stages
  all_info$state_specificity <- apply(cell_state_means[, all_info$TF, drop = FALSE], 2L, stats::sd)
  all_info$stage_range <- apply(unit_stage_matrix[all_info$TF, , drop = FALSE], 1L,
                                function(v) max(v) - min(v))
  n_candidates <- nrow(all_info)
  # Equal weight to state discrimination and composition-adjusted stage
  # variation; descriptive ranking, not significant differential regulation.
  all_info$state_rank_pct <- rank(all_info$state_specificity, ties.method = "average") / n_candidates
  all_info$stage_rank_pct <- rank(all_info$stage_range, ties.method = "average") / n_candidates
  all_info$priority_score <- (all_info$state_rank_pct + all_info$stage_rank_pct) / 2
  all_info <- all_info[order(-all_info$priority_score,
                             -all_info$stage_range, -all_info$state_specificity,
                             all_info$TF), , drop = FALSE]
  all_info$rank <- seq_len(nrow(all_info))
  utils::write.csv(all_info, file.path(OUT_DIR, "Fig5F_All_Candidate_TF_Ranking.csv"), row.names = FALSE)
  network_info <- utils::head(all_info, 50L)
  network_tf <- network_info$TF
  target_sets <- setNames(lapply(network_info$regulon,
                                  function(r) as_genes(regulon_list[[r]])), network_tf)
  pairs <- utils::combn(network_tf, 2L, simplify = FALSE)
  edges_all <- do.call(rbind, lapply(pairs, function(two) {
    shared <- length(intersect(target_sets[[two[1L]]], target_sets[[two[2L]]]))
    total <- length(union(target_sets[[two[1L]]], target_sets[[two[2L]]]))
    data.frame(from = two[1L], to = two[2L],
               n_shared = shared, jaccard = if (total) shared / total else 0)
  }))
  edges_all <- edges_all[order(-edges_all$jaccard, -edges_all$n_shared), , drop = FALSE]
  edges <- edges_all[edges_all$n_shared >= 5L & edges_all$jaccard >= 0.02, , drop = FALSE]
  if (!nrow(edges)) edges <- edges_all[edges_all$n_shared >= 1L, , drop = FALSE]
  edges <- utils::head(edges, 85L)  # selected strongest pairs, all 50 nodes shown

  # Stage colours come from independently averaged, state-centered AUCell.
  candidate_stage_means <- unit_stage_matrix[network_tf, stages, drop = FALSE]
  candidate_stage_range <- apply(candidate_stage_means, 1L, function(v) max(v) - min(v))
  node_stage <- stages[max.col(candidate_stage_means, ties.method = "first")]
  # Do not assign a visually decisive stage for essentially flat activity.
  node_stage[candidate_stage_range < 0.25] <- "No clear peak"
  nodes <- data.frame(
    name = network_tf,
    highest_relative_stage = factor(node_stage, levels = c(stages, "No clear peak")),
    n_targets = lengths(target_sets),
    state_specificity = network_info$state_specificity,
    stage_range = candidate_stage_range,
    rank = network_info$rank,
    regulon = network_info$regulon,
    stringsAsFactors = FALSE
  )
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  set.seed(20261023)
  layout <- if (igraph::ecount(graph) > 0L) {
    igraph::layout_with_fr(graph, weights = igraph::E(graph)$jaccard + 0.01)
  } else {
    igraph::layout_in_circle(graph)
  }
  node_positions <- data.frame(name = igraph::V(graph)$name,
                               x = layout[, 1L], y = layout[, 2L])
  nodes <- merge(nodes, node_positions, by = "name", sort = FALSE)
  if (igraph::ecount(graph) > 0L) {
    ends <- igraph::ends(graph, igraph::E(graph), names = TRUE)
    line_data <- data.frame(
      from = ends[, 1L], to = ends[, 2L],
      n_shared = igraph::E(graph)$n_shared,
      jaccard = igraph::E(graph)$jaccard,
      stringsAsFactors = FALSE
    )
    line_data$x <- node_positions$x[match(line_data$from, node_positions$name)]
    line_data$y <- node_positions$y[match(line_data$from, node_positions$name)]
    line_data$xend <- node_positions$x[match(line_data$to, node_positions$name)]
    line_data$yend <- node_positions$y[match(line_data$to, node_positions$name)]
  } else {
    line_data <- data.frame(x = numeric(), y = numeric(), xend = numeric(),
                            yend = numeric(), from = character(), to = character(),
                            n_shared = integer(), jaccard = numeric())
  }
  # Label informative hubs rather than drawing 50 overlapping strings.
  hub_strength <- setNames(rep(0, length(network_tf)), network_tf)
  if (nrow(edges)) {
    for (i in seq_len(nrow(edges))) {
      hub_strength[edges$from[i]] <- hub_strength[edges$from[i]] + edges$jaccard[i]
      hub_strength[edges$to[i]] <- hub_strength[edges$to[i]] + edges$jaccard[i]
    }
  }
  nodes$network_strength <- unname(hub_strength[nodes$name])
  label_tf <- utils::head(nodes$name[order(-nodes$network_strength, nodes$rank)], 16L)
  nodes$label <- ifelse(nodes$name %in% label_tf, nodes$name, "")
  pF <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = line_data,
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                                       linewidth = jaccard),
                          colour = "#A8B3B9", alpha = 0.46) +
    ggplot2::geom_point(data = nodes,
                        ggplot2::aes(x = x, y = y, fill = highest_relative_stage,
                                     size = state_specificity),
                        shape = 21, colour = "white", stroke = 0.65)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    pF <- pF + ggrepel::geom_text_repel(
      data = nodes[nzchar(nodes$label), , drop = FALSE],
      ggplot2::aes(x = x, y = y, label = label),
      seed = 20261023, max.overlaps = Inf, min.segment.length = 0,
      fontface = "bold", size = 3.7, colour = "#263742", box.padding = 0.4
    )
  } else {
    pF <- pF + ggplot2::geom_text(
      data = nodes[nzchar(nodes$label), , drop = FALSE],
      ggplot2::aes(x = x, y = y, label = label),
      vjust = -1.2, fontface = "bold", size = 3.2, colour = "#263742"
    )
  }
  pF <- pF +
    ggplot2::scale_fill_manual(values = c(palette_stage, "No clear peak" = "#ACB5BA"),
                               name = "Highest relative\nstage activity", drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(3.6, 8.5), name = "State specificity") +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 1.75), name = "Target overlap\n(Jaccard)") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = "F  |  Fifty candidate TF regulons and shared targets",
                  subtitle = "Ranked jointly by myeloid-state specificity and stage-associated variation",
                  caption = "50 TFs shown; only 16 hubs labeled. Edges indicate shared inferred targets, NOT proven TF-to-TF regulation.\nStage labels are descriptive and hypothesis-generating.") +
    ggplot2::theme_void(base_size = 11, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#58656D"),
      plot.caption = ggplot2::element_text(size = 8.5, colour = "#62717A", hjust = 0),
      legend.position = "right",
      plot.margin = ggplot2::margin(20, 40, 25, 20)
    )
  utils::write.csv(edges_all, file.path(OUT_DIR, "Fig5F_All_TF_Target_Overlap.csv"), row.names = FALSE)
  utils::write.csv(line_data, file.path(OUT_DIR, "Fig5F_Plotted_TF_Target_Network_Edges.csv"), row.names = FALSE)
  utils::write.csv(nodes, file.path(OUT_DIR, "Fig5F_TF_Target_Network_Nodes.csv"), row.names = FALSE)
  save_fig(pF, "Fig5F_Top50_SCENIC_TF_Target_Network", 13, 9.5)

  cat("\n=====================================================\n")
  cat("FIG5 A-F GENERATED FROM 1529 DECONTX-CORRECTED CELLS\n")
  cat("Cached AUCell reused; no SCENIC recalculation.\n")
  cat("No significance claims from individual cells.\n")
  cat("Output directory: ", OUT_DIR, "\n", sep = "")
  cat("=====================================================\n")
  print(pA)
  print(pB)
  print(pC)
  print(pD)
  print(pE)
  print(pF)
})



}

# END OF FILE.
# Before E3 release: replace all GAPs with verified historical scripts/inputs,
# review every legacy path, run a clean-session reproducibility check,
# then distribute only the tested version + real dependencies and input provenance.
