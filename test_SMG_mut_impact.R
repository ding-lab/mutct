# Yige Wu @ WashU 2018 Aug
## Revised by Yize Li, Song Cao
## test mutation impact on protein/phosphorylation within kinase-substrate pairs or protein complex pairs
## using wilcox test
## reference: https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/

# source ------------------------------------------------------------------
setwd(dir = "/Users/scao/Desktop/all/projects/cptac3/mutation_effect/scripts/R/Mutational_Impact/")

source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sample_type <- "tumor"
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "num", "meddiff_bottom", "meddiff_upper", "meddiff")
num_control_thres <- 4
num_genoalt_thres <- 4

# cancer types to process -------------------------------------------------
cancers2process <- c("LSCC")

# expression data to process ----------------------------------------------
affected_exp_types2process <- c("PRO")

# variant types to test ---------------------------------------------------
## include all nonsynonymous mutations or just missense/truncation
# variant_classes2process <- c("fusion")

variant_classes2process <- c("truncation", "missense", "not_silent")

# input protein pair table ------------------------------------------------
pair_tab_annotated <- fread(input = "./protein_pair_table.txt", data.table = F)

#SMG_list <- c("DNAH14", "DNAH9", "FMN2", "METTL4", "MYT1L", "OBSCN", "PLXNA2", "RYR3", "TAF1L", "MUC16", "RNF43", "RYR2", "ZFHX4", "CTNNA2", "MUC4", "TTN", "CDKN2A", "SMAD4", "TP53", "KRAS", "ARID1A", "TGFBR2", "GNAS", "RREB1", "PBRM1", "BRAF", "CTNNB1", "U2AF1", "KDM6A")

gene_list <- "~/Desktop/all/projects/data_source/smg.bailey.cell.lusc.tsv"
gene_list <- read.table(gene_list)

SMG_list <- gene_list$V1

pair_tab <- pair_tab_annotated[pair_tab_annotated$GENE %in% SMG_list, c("GENE", "SUB_GENE")]

# Business  ---------------------------------------------------------------
# affected_exp_type <- "Pro"
# variant_class <- "not_silent"
# cancer <- "LSCC"
for (affected_exp_type in affected_exp_types2process) {
  for (variant_class in variant_classes2process) {
   for (cancer in cancers2process) {
      fn <-  paste0(cancer, "_", variant_class, "_mut_impact_", affected_exp_type , "_tab.v3.2.txt")
      if (!file.exists(fn) || file.exists(fn)) {
        if (affected_exp_type == "PRO") {
          affected_exp_data <- fread("/Users/scao/Desktop/all/projects/cptac3/analysis/LSCC/data/v3.2p/prot.lscc.allgene.tumor.tsv", data.table = F)
          affected_exp_head <- data.frame(SUBSTRATE = affected_exp_data$Gene)
        }

        ## input mutation matrix
        maf <- fread(input = "/Users/scao/Desktop/all/projects/cptac3/analysis/LSCC/data/v3.2p/somatic-lscc-v3.2.maf", data.table = F)
        mut_mat <- generate_somatic_mutation_matrix(pair_tab = pair_tab$GENE, maf = maf)
        maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]

        ## get the overlap of the participant IDs
        partIDs <- colnames(affected_exp_data)
        length(partIDs)
        length(colnames(mut_mat))
        partIDs_overlap <- intersect(partIDs, colnames(mut_mat))
        length(partIDs_overlap)

        ## the number of mutations
        if (variant_class == "truncation") {
          mut_count <- sapply(1:nrow(mut_mat), FUN = function(i, mut_mat) sum(grepl(x = mut_mat[i,], pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del|Nonstop_Mutation|Splice_Site")), mut_mat[,-1])
        }
        if (variant_class == "missense") {
          mut_count <- sapply(1:nrow(mut_mat), FUN = function(i, mut_mat) sum(grepl(x = mut_mat[i,], pattern = "Missense_Mutation|In_Frame_Ins|In_Frame_Del")), mut_mat[,-1])
        }
        if (variant_class == "not_silent") {
          mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent")))
        }
        if (variant_class == "fusion") {
          mut_count <- sapply(1:nrow(mut_mat), FUN = function(i, mut_mat) sum(grepl(x = mut_mat[i,], pattern = "Fusion")), mut_mat[,-1])
        }
        names(mut_count) <- mut_mat$Hugo_Symbol
        # print(mut_count)
        # readline(prompt="Press [enter] to continue")
        ## get the genes to test
        genes_num_genoalt_thresholded <- names(mut_count)[mut_count >= num_genoalt_thres]

        # print(length(genes_num_genoalt_thresholded))
        # readline(prompt="Press [enter] to continue")

        ## initiate identifier columns
        pair_tab_genoalt_thresholded <- pair_tab[pair_tab$GENE %in% genes_num_genoalt_thresholded,]
        pair_tab_genoalt_thresholded
        # print(nrow(pair_tab_genoalt_thresholded))

        if (affected_exp_type %in% c("PRO", "RNA")) {
          pair_tab_sites <- pair_tab_genoalt_thresholded
          pair_tab_sites$SUB_MOD_RSD <- affected_exp_type
        } else {
          pair_tab_sites <- merge(pair_tab_genoalt_thresholded, affected_exp_head[, c("SUBSTRATE", "SUB_MOD_RSD")], by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"))
        }
        pair_tab_sites <- pair_tab_sites[!is.na(pair_tab_sites$GENE) & !is.na(pair_tab_sites$SUB_MOD_RSD),]
        pair_tab_sites <- pair_tab_sites[order(pair_tab_sites$GENE, pair_tab_sites$SUB_GENE),]
        print(nrow(pair_tab_sites))

        # readline(prompt="Press [enter] to continue")
        ## initiate value columns
        num_control <- vector(mode = "numeric", length = nrow(pair_tab_sites))
        num_mut <- num_control; meddiff_mut <- num_control; p_mut <- num_control; meddiff_bottom_mut <- num_control; meddiff_upper_mut <- num_control


        for (enzyme in unique(pair_tab_sites$GENE)) {
          mut_mat_en <- mut_mat[mut_mat$Hugo_Symbol == enzyme,]
          if (nrow(mut_mat_en) > 0){
            if (variant_class == "not_silent") {
              mut_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] != "") & (mut_mat_en[, partIDs_overlap] != "Silent")]
            }
            if (variant_class == "missense") {
              mut_partIDs <- partIDs_overlap[grepl(x = mut_mat_en[, partIDs_overlap], pattern = "Missense_Mutation|In_Frame_Ins|In_Frame_Del")]
            }
            if (variant_class == "truncation") {
              mut_partIDs <- partIDs_overlap[grepl(x = mut_mat_en[, partIDs_overlap], pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del|Nonstop_Mutation|Splice_Site")]
            }
            if (variant_class == "fusion") {
              mut_partIDs <- partIDs_overlap[grepl(x = mut_mat_en[, partIDs_overlap], pattern = "Fusion")]
            }
            

          } else {
            mut_partIDs <- NULL
          }


          mut_cnv_partIDs <- mut_partIDs
          mut_cnv_partIDs
          control_partIDs <- partIDs_overlap[(mut_mat_en[,partIDs_overlap] == "") | (mut_mat_en[, partIDs_overlap] == "Silent")]

          # readline(prompt="Press [enter] to continue")
          print(paste0("enzyme:", enzyme))
          for (substrate in unique(pair_tab_sites$SUB_GENE[pair_tab_sites$GENE == enzyme])) {
            print(paste0("substrate:", substrate))

            for (site in unique(pair_tab_sites$SUB_MOD_RSD[pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate])) {
              print(paste0("site:", site))

              if (affected_exp_type  %in% c("PRO", "RNA")) {
                affected_exp_altered <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate,]
                i <- (pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate & pair_tab_sites$SUB_MOD_RSD == site)
              } else {
                affected_exp_altered <- affected_exp_data[affected_exp_head$SUBSTRATE == substrate & affected_exp_head$SUB_MOD_RSD == site,]
                i <- (pair_tab_sites$GENE == enzyme & pair_tab_sites$SUB_GENE == substrate & pair_tab_sites$SUB_MOD_RSD == site)

              }
              
              affected_exp_altered <- affected_exp_altered[1,]
              affected_exp_control <- affected_exp_altered[, intersect(control_partIDs, colnames(affected_exp_altered))]
              affected_exp_control <- affected_exp_control[!is.na(affected_exp_control)]

              num_control[i] <- length(affected_exp_control)
              # print(affected_exp_control)
              if (num_control[i] >= num_control_thres) {
                mut_pho <- affected_exp_altered[, mut_partIDs]
                mut_pho <- mut_pho[!is.na(mut_pho)]
                num_mut[i] <- length(mut_pho)
                if (num_mut[i] > 0) {
                  meddiff_mut[i] <- median(mut_pho) - median(affected_exp_control)
                  if (num_mut[i] >= num_genoalt_thres){
                    # readline(prompt="Press [enter] to continue: mut")
                    # print(mut_pho)
                    # readline(prompt="Press [enter] to continue: control")
                    # print(affected_exp_control)
                    # 
                    stat <- wilcox.test(x = mut_pho, y = affected_exp_control, conf.int = T)
                    p_mut[i] <- stat$p.value
                    meddiff_bottom_mut[i] <- stat$conf.int[1]
                    meddiff_upper_mut[i] <- stat$conf.int[2]
                  }
                }
              }
              # print(which(i))
            }
          }

        }
        # print(p_mut)
        mut_cnv_tab <- pair_tab_sites
        mut_cnv_tab <- cbind(mut_cnv_tab,
                             data.frame(num = num_mut,
                                        num_control = num_control,
                                        p = p_mut,
                                        p_fdr = p.adjust(p_mut, method = "fdr", n = length(p_mut)),
                                        p_bonferroni = p.adjust(p_mut, method = "bonferroni", n = length(p_mut)),
                                        meddiff_bottom = meddiff_bottom_mut,
                                        meddiff_upper = meddiff_upper_mut,
                                        meddiff = meddiff_mut))
        # print(p_mut)
        # readline(prompt="Press [enter] to continue: p_mut")
        mut_cnv_tab$cancer <- cancer
        mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
        ## annotate cis and trans
        mut_cnv_tab$SELF <- "trans"
        mut_cnv_tab$SELF[as.vector(mut_cnv_tab$GENE) == as.vector(mut_cnv_tab$SUB_GENE)] <- "cis"
        mut_cnv_tab$genoalt_type <- "mut"

        ## clean up

        write.table(x = mut_cnv_tab, file = fn, col.names = T, row.names = F, quote = F, sep = "\t")

      }
    }
  }
}
