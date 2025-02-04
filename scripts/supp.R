# R functions

nicefy_gene_results <- function(input, genelist) {
  temp_table <- input
  genes_in_input <- subset(genelist,
                           (genelist$Gene_stable_ID %in% input$Gene_stable_ID))
  
  out_table <- merge(
    temp_table,
    genes_in_input,
    by.x = c("Gene_stable_ID"),
    by.y = c("Gene_stable_ID"),
    all = TRUE
  ) %>%
    relocate(Gene_name, .after = Gene_stable_ID)
  
  return(out_table)
}

intersectionalize <- function(set1, set2, set3) {
  # Selects the intersecting Gene_stable_ID's and extracts log2Foldchange from each experient for later merging.
  tmp1 <- set1 %>% subset(Gene_stable_ID %in% set2$Gene_stable_ID) %>% subset(Gene_stable_ID %in% set3$Gene_stable_ID) %>% dplyr::select(Gene_stable_ID, Gene_name, log2FoldChange)
  tmp2 <- set2 %>% subset(Gene_stable_ID %in% tmp1$Gene_stable_ID) %>% dplyr::select(Gene_stable_ID, log2FoldChange)
  tmp3 <- set3 %>% subset(Gene_stable_ID %in% tmp1$Gene_stable_ID) %>% dplyr::select(Gene_stable_ID, log2FoldChange)
  
  # Merges the three experiments and sets colnames
  tmp4 <- merge(tmp1, tmp2, by = "Gene_stable_ID")
  tmp4 <- merge(tmp4, tmp3, by = "Gene_stable_ID")
  colnames(tmp4) <- c("Gene_stable_ID", "Gene_name", "l2FC_1", "l2FC_2", "l2FC_3")
  
  # Calculates how many of the intersecting genes have the same direction (signs).
  # Also calculates the mean and std
  tmp4 <- tmp4 %>% mutate(signs = rowSums(sign(tmp4[, 3:5])))
  tmp4 <- tmp4 %>% mutate(m = if_else(abs(signs) == 3, apply(tmp4[, 3:5], 1, mean), NA))
  tmp4 <- tmp4 %>% mutate(s = if_else(abs(signs) == 3, apply(tmp4[, 3:5], 1, sd), NA))
  
  return(tmp4)
}

FGSEA_wrap <- function(input, gscGO) {
  tmp <- input %>% dplyr::filter(!is.na(padj) &
                                   !is.na(pvalue) & baseMean > 0)
  Pval <- tmp$pvalue
  Pval[Pval == 0] <- min(Pval[Pval != 0])
  names(Pval) <- tmp$Gene_stable_ID
  FC <- tmp$log2FoldChange
  names(FC) <- tmp$Gene_stable_ID
  gsaRes <- runGSA(
    -log10(Pval) * sign(FC),
    geneSetStat = "fgsea",
    gsc = gscGO,
    gsSizeLim = c(10, 200),
    nPerm = 10000
  )
  return(gsaRes)
}

gene_checker_deluxe <- function(input1,
                                input2,
                                input3,
                                gene,
                                contrast,
                                lims = F) {
  if (gene %in% input1$Gene_name) {
    tmp1 <- input1 %>% subset(Gene_name == gene)
    tmp2 <- input2 %>% subset(Gene_name == gene)
    tmp3 <- input3 %>% subset(Gene_name == gene)
  } else {
    tmp1 <- input1 %>% subset(Gene_stable_ID == gene)
    tmp2 <- input2 %>% subset(Gene_stable_ID == gene)
    tmp3 <- input3 %>% subset(Gene_stable_ID == gene)
  }
  
  contr <- case_when(
    contrast == "log" ~ "Pre",
    contrast == "PDS" ~ "Post",
    contrast == "Ref" ~ "^R",
    contrast == "TDA" ~ "^T"
  )
  
  tmp1 <- tmp1[, grepl(contr, colnames(tmp1))]
  tmp2 <- tmp2[, grepl(contr, colnames(tmp2))]
  tmp3 <- tmp3[, grepl(contr, colnames(tmp3))]
  
  rownames(tmp1) <- "Gene"
  rownames(tmp2) <- "Gene"
  rownames(tmp3) <- "Gene"
  
  colnames(tmp1) <- paste0(colnames(tmp1), "1")
  colnames(tmp2) <- paste0(colnames(tmp2), "2")
  colnames(tmp3) <- paste0(colnames(tmp3), "3")
  
  slist1 <- colnames(tmp1)
  slist2 <- colnames(tmp2)
  slist3 <- colnames(tmp3)
  
  tmp4 <- gsub("[0-9]", ".", colnames(tmp1)) %>% paste0(., "1") %>% strsplit("[.]") %>% as.data.frame()
  tmp5 <- gsub("[0-9]", ".", colnames(tmp2)) %>% paste0(., "2") %>% strsplit("[.]") %>% as.data.frame()
  tmp6 <- gsub("[0-9]", ".", colnames(tmp3)) %>% paste0(., "3") %>% strsplit("[.]") %>% as.data.frame()
  
  colnames(tmp4) <- slist1
  colnames(tmp5) <- slist2
  colnames(tmp6) <- slist3
  
  rownames(tmp4) <- c("Strain", "Phase" , "Experiment")
  rownames(tmp5) <- c("Strain", "Phase" , "Experiment")
  rownames(tmp6) <- c("Strain", "Phase" , "Experiment")
  
  tmp1 <- rbind(tmp1, tmp4) %>% t()
  tmp2 <- rbind(tmp2, tmp5) %>% t()
  tmp3 <- rbind(tmp3, tmp6) %>% t()
  
  tmp1 <- tmp1 %>% as.data.frame()
  tmp2 <- tmp2 %>% as.data.frame()
  tmp3 <- tmp3 %>% as.data.frame()
  
  tmp <- rbind(tmp1, tmp2, tmp3)
  
  tmp$Strain <- factor(tmp$Strain, levels = c("R", "T"), ordered = T)
  tmp$Phase <- factor(tmp$Phase,
                      levels = c("Pre", "Post"),
                      ordered = T)
  tmp$Experiment <- factor(tmp$Experiment,
                           levels = c("1", "2", "3"),
                           ordered = T)
  tmp$Gene <- as.numeric(tmp$Gene)
  
  if (contrast %in% c("log", "PDS")) {
    p <- ggplot(tmp, aes(x = Experiment, y = Gene, fill = Strain))
    groupcol <- c("#8ECA4F", "#E20E6F")
    grouplab <- c("Reference", "tda1∆")
  } else {
    p <- ggplot(tmp, aes(x = Experiment, y = Gene, fill = Phase))
    groupcol <- c("#E3D455", "#00588E")
    grouplab <- c("log-phase", "PDS-phase")
  }
  
  p <- p + geom_boxplot() + theme_classic() + ggtitle(gene) + scale_fill_manual(values = groupcol, labels = grouplab) + {
    if (any(lims != F))
      ylim(lims[1], lims[2])
  } +
    labs(y = "Normalized gene count") +
    # guides(fill = )
    theme(
      plot.margin = unit(c(1.5, 1, 1.5, 1), "lines"),
      title = element_text(size = 20),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.key.size = unit(2, "lines")
    )
  return(p)
}

gsdot <- function(ref1,
                  ref2,
                  ref3,
                  mut1,
                  mut2,
                  mut3,
                  exps = c("1", "2", "3"),
                  strs = c("Reference", "tda1∆"),
                  plim = 0.01) {
  TDAup1 <- data.frame(
    GO = names(mut1$gsc),
    Gene_ratio = unname(mut1$nGenesUp) / unname(mut1$nGenesTot),
    padj = mut1$pAdjDistinctDirUp,
    Experiment = exps[1],
    Strain = strs[2],
    Direction = "Up"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  TDAup2 <- data.frame(
    GO = names(mut2$gsc),
    Gene_ratio = unname(mut2$nGenesUp) / unname(mut2$nGenesTot),
    padj = mut2$pAdjDistinctDirUp,
    Experiment = exps[2],
    Strain = strs[2],
    Direction = "Up"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  TDAup3 <- data.frame(
    GO = names(mut3$gsc),
    Gene_ratio = unname(mut3$nGenesUp) / unname(mut3$nGenesTot),
    padj = mut3$pAdjDistinctDirUp,
    Experiment = exps[3],
    Strain = strs[2],
    Direction = "Up"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  Refup1 <- data.frame(
    GO = names(ref1$gsc),
    Gene_ratio = unname(ref1$nGenesUp) / unname(ref1$nGenesTot),
    padj = ref1$pAdjDistinctDirUp,
    Experiment = exps[1],
    Strain = strs[1],
    Direction = "Up"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  Refup2 <- data.frame(
    GO = names(ref2$gsc),
    Gene_ratio = unname(ref2$nGenesUp) / unname(ref2$nGenesTot),
    padj = ref2$pAdjDistinctDirUp,
    Experiment = exps[2],
    Strain = strs[1],
    Direction = "Up"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  Refup3 <- data.frame(
    GO = names(ref3$gsc),
    Gene_ratio = unname(ref3$nGenesUp) / unname(ref3$nGenesTot),
    padj = ref3$pAdjDistinctDirUp,
    Experiment = exps[3],
    Strain = strs[1],
    Direction = "Up"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  TDAdn1 <- data.frame(
    GO = names(mut1$gsc),
    Gene_ratio = unname(mut1$nGenesDn) / unname(mut1$nGenesTot),
    padj = mut1$pAdjDistinctDirDn,
    Experiment = exps[1],
    Strain = strs[2],
    Direction = "Down"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  TDAdn2 <- data.frame(
    GO = names(mut2$gsc),
    Gene_ratio = unname(mut2$nGenesDn) / unname(mut2$nGenesTot),
    padj = mut2$pAdjDistinctDirDn,
    Experiment = exps[2],
    Strain = strs[2],
    Direction = "Down"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  TDAdn3 <- data.frame(
    GO = names(mut3$gsc),
    Gene_ratio = unname(mut3$nGenesDn) / unname(mut3$nGenesTot),
    padj = mut3$pAdjDistinctDirDn,
    Experiment = exps[3],
    Strain = strs[2],
    Direction = "Down"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  Refdn1 <- data.frame(
    GO = names(ref1$gsc),
    Gene_ratio = unname(ref1$nGenesDn) / unname(ref1$nGenesTot),
    padj = ref1$pAdjDistinctDirDn,
    Experiment = exps[1],
    Strain = strs[1],
    Direction = "Down"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  Refdn2 <- data.frame(
    GO = names(ref2$gsc),
    Gene_ratio = unname(ref2$nGenesDn) / unname(ref2$nGenesTot),
    padj = ref2$pAdjDistinctDirDn,
    Experiment = exps[2],
    Strain = strs[1],
    Direction = "Down"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  Refdn3 <- data.frame(
    GO = names(ref3$gsc),
    Gene_ratio = unname(ref3$nGenesDn) / unname(ref3$nGenesTot),
    padj = ref3$pAdjDistinctDirDn,
    Experiment = exps[3],
    Strain = strs[1],
    Direction = "Down"
  ) %>% na.omit() %>% dplyr::filter(padj < plim)
  
  allofit <- rbind(
    TDAup1,
    TDAup2,
    TDAup3,
    Refup1,
    Refup2,
    Refup3,
    TDAdn1,
    TDAdn2,
    TDAdn3,
    Refdn1,
    Refdn2,
    Refdn3
  )
  
  allofit$Direction <- factor(allofit$Direction, levels = c("Up", "Down"))
  
  p <- ggplot(allofit,
              aes(
                x = Experiment,
                y = GO,
                size = Gene_ratio,
                fill = -log10(padj)
              )) +
    geom_point(pch = 21) +
    scale_y_discrete(position = "right") +
    facet_grid(
      Direction ~ Strain,
      scales = "free",
      space = "free",
      switch = "y"
    ) +
    scale_fill_gradient2(
      mid = "white",
      high = "purple",
      midpoint = 2,
      limits = c(2, 3)
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 15, color = "black"),
      axis.y.text = element_text(color = "black"),
      legend.position = "left"
    ) +
    labs(y = "GO-terms", fill = "-log10(FDR)", size = "Gene Ratio")
  
  return(p)
  
}


customheat <- function(counts, experiment, samples, set, legend = F) {
  tmp <- matrix(, nrow = 0, ncol = dim(counts)[2])
  for (i in unique(set$GSname)) {
    ttmp <- set %>% subset(GSname == i) %>% pull(Ensembl)
    if (length(ttmp) >= 5) {
      tttmp <- counts[rownames(counts) %in% ttmp, ]
      tmp <- rbind(tmp, tttmp)
    }
  }
  tmp <- tmp %>% na.omit() %>% apply(1, scale) %>% t() %>% na.omit()
  colnames(tmp) <- samples$Sample
  anno <- HeatmapAnnotation(
    Experiment = rep(experiment, nrow(samples)),
    Strain = samples %>% arrange(Phase) %>% arrange(Strain) %>% pull(Strain),
    Phase = samples %>% arrange(Phase) %>% arrange(Strain) %>% pull(Phase),
    col = list(
      Experiment = c(
        "1" = "#F58868",
        "2" = "#B4FFF5",
        "3" = "#A154FF"
      ),
      Strain = c(Ref = "#8ECA4F", "tda1∆" = "#E20E6F"),
      Phase = c(log = "#E3D455", PDS = "#227EB7")
    ),
    border = T,
    gp = gpar(col = "black"),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 15, fontface = "bold"),
      labels_gp = gpar(fontsize = 12)
    )
  )
  
  h <- Heatmap(
    tmp,
    cluster_rows = T,
    cluster_columns = F,
    top_annotation = anno,
    show_row_names = F,
    show_heatmap_legend = legend,
    column_names_gp = gpar(fontsize = 10),
    col = colorRamp2(c(-2, 0, 2), c(
      plasma(100)[1], plasma(100)[50], plasma(100)[100]
    )),
    row_names_gp = grid::gpar(fontsize = 4),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 15, fontface = "bold"),
      labels_gp = gpar(fontsize = 12),
      title = "Z-score",
      at = seq(-2, 2, 1)
    ),
    width = unit(ncol(counts) * 0.3, "cm")
  )
  return(h)
}

mycomb <- function(input) {
  tmp <- list_to_matrix(input)
  tmp <- data.frame(tmp)
  tmp <- tmp %>% mutate(
    ol = l * (p + 1) %% 2 * (r + 1) %% 2 * (t + 1) %% 2,
    op = (l + 1) %% 2 * p * (r + 1) %% 2 * (t + 1) %% 2,
    or = (l + 1) %% 2 * (p + 1) %% 2 * r * (t + 1) %% 2,
    ot = (l + 1) %% 2 * (p + 1) %% 2 * (r + 1) %% 2 * t,
    lp = l * p * (r + 1) %% 2 * (t + 1) %% 2,
    lr = l * (p + 1) %% 2 * r * (t + 1) %% 2,
    lt = l * (p + 1) %% 2 * (r + 1) %% 2 * t,
    pr = (l + 1) %% 2 * p * r * (t + 1) %% 2,
    pt = (l + 1) %% 2 * p * (r + 1) %% 2 * t,
    rt = (l + 1) %% 2 * (p + 1) %% 2 * r * t,
    lpr = l * p * r * (t + 1) %% 2,
    lpt = l * p * (r + 1) %% 2 * t,
    lrt = l * (p + 1) %% 2 * r * t,
    prt = (l + 1) %% 2 * p * r * t,
    lprt = l * p * r * t
  )
  return(tmp)
}

btest <- function (x1,
                   n1,
                   x2,
                   n2,
                   alternative = c("two.sided", "less", "greater"),
                   or = NULL,
                   conf.int = FALSE,
                   conf.level = 0.95,
                   midp = FALSE,
                   tsmethod = c("central", "minlike"),
                   control = ucControl())
{
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  if (midp) {
    T <- function(X1, N1, X2, N2, delta0 = 1) {
      m <- length(X1)
      pval <- rep(NA, m)
      for (i in 1:m) {
        pval[i] <- exact2x2(
          matrix(c(X2[i], N2 - X2[i], X1[i], N1 - X1[i]), 2, 2),
          or = delta0,
          conf.int = FALSE,
          midp = TRUE,
          alternative = "less"
        )$p.value
      }
      pval
    }
  }
  else if (alternative == "less") {
    T <- function(X1, N1, X2, N2, delta0 = 1) {
      m <- length(X1)
      pval <- rep(NA, m)
      for (i in 1:m) {
        pval[i] <- exact2x2(
          matrix(c(X2[i], N2 - X2[i], X1[i], N1 - X1[i]), 2, 2),
          or = delta0,
          conf.int = FALSE,
          midp = FALSE,
          alternative = "less"
        )$p.value
      }
      pval
    }
  }
  else if (alternative == "greater") {
    T <- function(X1, N1, X2, N2, delta0 = 1) {
      m <- length(X1)
      pval <- rep(NA, m)
      for (i in 1:m) {
        pval[i] <- exact2x2(
          matrix(c(X2[i], N2 - X2[i], X1[i], N1 - X1[i]), 2, 2),
          or = delta0,
          conf.int = FALSE,
          midp = FALSE,
          alternative = "greater"
        )$p.value
      }
      1 - pval
    }
  }
  else if (alternative == "two.sided" & tsmethod == "minlike") {
    T <- function(X1, N1, X2, N2, delta0 = 1) {
      m <- length(X1)
      pval <- rep(NA, m)
      for (i in 1:m) {
        cat("\r",
            paste0(
              "Minlike: ",
              i,
              " out of ",
              m,
              " completed ",
              "(",
              round(i / m * 100, 0),
              "%)."
            ))
        pval[i] <- fisher.test(
          matrix(c(X2[i], N2 -
                     X2[i], X1[i], N1 - X1[i]), 2, 2),
          or = delta0,
          conf.int = FALSE,
          alternative = "two.sided"
        )$p.value
      }
      1 - pval
    }
    out <- uncondExact2x2(
      x1,
      n1,
      x2,
      n2,
      alternative = alternative,
      nullparm = or,
      parmtype = "oddsratio",
      conf.int = conf.int,
      conf.level = conf.level,
      midp = midp,
      method = "user",
      tsmethod = "square",
      control = control,
      Tfunc = T
    )
    out$method <- "Boschloo's test"
    return(out)
  }
  else if (alternative == "two.sided" & tsmethod == "central") {
    conf.level <- 1 - (1 - conf.level) / 2
    T <- function(X1, N1, X2, N2, delta0 = 1) {
      m <- length(X1)
      pval <- rep(NA, m)
      for (i in 1:m) {
        cat("\r", paste0("Less:", i, " out of ", m, " completed."))
        pval[i] <- exact2x2(
          matrix(c(X2[i], N2 - X2[i], X1[i], N1 - X1[i]), 2, 2),
          or = delta0,
          conf.int = FALSE,
          midp = FALSE,
          alternative = "less"
        )$p.value
      }
      pval
    }
    outLess <- uncondExact2x2(
      x1,
      n1,
      x2,
      n2,
      alternative = "less",
      nullparm = or,
      parmtype = "oddsratio",
      conf.int = conf.int,
      conf.level = conf.level,
      midp = midp,
      method = "user",
      control = control,
      Tfunc = T
    )
    T <- function(X1, N1, X2, N2, delta0 = 1) {
      m <- length(X1)
      pval <- rep(NA, m)
      for (i in 1:m) {
        cat("\r", paste0("Greater:", i, " out of ", m, " completed."))
        pval[i] <- exact2x2(
          matrix(c(X2[i], N2 - X2[i], X1[i], N1 - X1[i]), 2, 2),
          or = delta0,
          conf.int = FALSE,
          midp = FALSE,
          alternative = "greater"
        )$p.value
      }
      1 - pval
    }
    outGreater <- uncondExact2x2(
      x1,
      n1,
      x2,
      n2,
      alternative = "greater",
      nullparm = or,
      parmtype = "oddsratio",
      conf.int = conf.int,
      conf.level = conf.level,
      midp = midp,
      method = "user",
      control = control,
      Tfunc = T
    )
    out <- outLess
    out$p.value <- min(c(1, 2 * outLess$p.value, 2 * outGreater$p.value))
    out$conf.int <- c(outGreater$conf.int[1], outLess$conf.int[2])
    out$method <- "Boschloo's test"
    return(out)
  }
  uncondMethod <- "user"
  out <- uncondExact2x2(
    x1,
    n1,
    x2,
    n2,
    alternative = alternative,
    nullparm = or,
    parmtype = "oddsratio",
    conf.int = conf.int,
    conf.level = conf.level,
    midp = midp,
    method = uncondMethod,
    control = control,
    Tfunc = T
  )
  out$method <- "Boschloo's test"
  if (tsmethod == "central" & alternative == "two.sided")
    out$method <- "Central Boschloo's test"
  out
}

megatest <- function(reslist, GS, seed) {
  set.seed(seed)
  gpes <- list()
  for (i in seq(reslist)) {
    g <- FGSEA_wrap(reslist[[i]], GS)
    gpes[[i]] <- g$pAdjDistinctDirUp
  }
  
  n <- nrow(reslist[[1]])
  out <- list()
  for (i in seq(GS$gsc)) {
    out[[i]] <- data.frame(
      Method = c(rep("Fisher", 3), rep("Boschloo", 3), rep("FGSEA", 3)),
      Experiment = rep(c(1, 2, 3), 3),
      log = rep(NA, 9),
      PDS = rep(NA, 9)
    )
    
    genelist <- GS$gsc[i]
    m <- length(genelist[[1]])
    for (j in seq(reslist)) {
      k <- reslist[[j]] %>% subset(padj < 0.05)
      q <- k %>% subset(Gene_stable_ID %in% genelist[[1]])
      
      k <- nrow(k)
      q <- nrow(q)
      
      print(paste0(
        "Currently doing ",
        names(GS$gsc)[i],
        ". Dataset ",
        j,
        " out of ",
        length(reslist),
        "."
      ))
      fish <- fisher.test(matrix(c(q, k - q, m - q, n - m - k + q), nr = 2), alternative = "g")$p.value
      bosch <- btest(k - q , n - m, q , m, tsmethod = "minlike")$p.value
      
      r <- j  %% 3 + 3 * !(j %% 3)
      c <-  floor((j - 1) / 3) + 3
      
      out[[i]][r + 3 * 0, c] <- fish
      out[[i]][r + 3 * 1, c] <- bosch
      out[[i]][r + 3 * 2, c] <- gpes[[j]][i]
    }
  }
  names(out) <- names(GS$gsc)
  return(out)
}



save_r_object <- function(object, path) {
  saveRDS(object, file = paste0(path, substitute(object), ".rds"))
}


xlsx_wrapper <- function(..., filename = "") {
  if (filename == "") {
    print("Please provide a filename")
    quit()
  }
  
  dfs <- list(...)
  wb = createWorkbook()
  i <- 1
  
  for (d in dfs) {
    sheet = createSheet(wb, paste0("Experiment ", i))
    addDataFrame(d, sheet = sheet, row.names = F)
    i <- i + 1
  }
  
  saveWorkbook(wb, filename)
}


export2svg <- function(object, size, filename = paste0(deparse(substitute(object)), ".svg")) {
  w <- size[1]
  h <- size[2]
  
  svglite(filename, width = w, height = h)
  try(plot(object), silent = T)
  dev.off()
}