# ============================================================
# ClinFold: Structural Impact Viewer for Genetic Variants
# Complete Shiny Application
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tibble)
  library(DT)
  library(bio3d)
  library(purrr)
  library(glue)
  library(stringr)
  library(ggplot2)
  library(plotly)
  library(NGLVieweR)
})

# ============================================================
# GLOBAL CONSTANTS
# ============================================================
UA <- "R httr (academic use) email: UO257655@uniovi.es"

# Complete amino acid mappings
base_map <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Glu="E", Gln="Q", Gly="G", 
  His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S", 
  Thr="T", Trp="W", Tyr="Y", Val="V"
)

AA_3TO1 <- c(
  setNames(base_map, tolower(names(base_map))),
  setNames(base_map, toupper(names(base_map))),
  base_map
)

AA_VALID <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", 
              "F", "P", "S", "T", "W", "Y", "V")

# Amino acid properties for scoring
AA_PROPERTIES <- list(
  # Hydrophobicity (Kyte-Doolittle scale, normalized 0-1)
  hydrophobicity = c(A=0.7, R=0.0, N=0.1, D=0.1, C=0.8, E=0.1, Q=0.1, G=0.5, 
                     H=0.2, I=1.0, L=0.9, K=0.0, M=0.7, F=0.9, P=0.3, S=0.3, 
                     T=0.4, W=0.6, Y=0.5, V=0.9),
  # Size (relative, 0-1)
  size = c(A=0.2, R=0.9, N=0.4, D=0.4, C=0.3, E=0.5, Q=0.5, G=0.0, H=0.6, 
           I=0.6, L=0.6, K=0.7, M=0.6, F=0.8, P=0.4, S=0.2, T=0.4, W=1.0, 
           Y=0.8, V=0.5),
  # Charge at pH 7
  charge = c(A=0, R=1, N=0, D=-1, C=0, E=-1, Q=0, G=0, H=0.1, I=0, L=0, K=1, 
             M=0, F=0, P=0, S=0, T=0, W=0, Y=0, V=0)
)

ALL_AA <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL")

`%||%` <- function(a, b) if (is.null(a)) b else a

DOMAIN_COLORS <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
                   "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62")

# ============================================================
# API FUNCTIONS
# ============================================================
safe_get <- function(url, timeout_secs = 30) {
  tryCatch({
    res <- httr::GET(url, httr::timeout(timeout_secs), httr::user_agent(UA))
    if (httr::http_error(res)) return(NULL)
    res
  }, error = function(e) NULL)
}

search_uniprot <- function(gene) {
  query <- paste0("gene:", gene, " AND organism_id:9606 AND reviewed:true")
  url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=", URLencode(query), "&format=json&size=5")
  
  r <- safe_get(url)
  if (is.null(r)) return(NULL)
  
  json <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  if (length(json$results) == 0) return(NULL)
  
  ent <- json$results[[1]]
  list(
    accession = ent$primaryAccession,
    id = ent$uniProtkbId,
    name = tryCatch(ent$proteinDescription$recommendedName$fullName$value, error = function(e) NA),
    gene = tryCatch(ent$genes[[1]]$geneName$value, error = function(e) gene),
    length = as.integer(ent$sequence$length),
    features = tryCatch({
      purrr::map_df(ent$features, function(f) {
        tibble(
          type = f$type %||% NA,
          description = f$description %||% NA,
          start = as.integer(f$location$start$value %||% NA),
          end = as.integer(f$location$end$value %||% NA)
        )
      })
    }, error = function(e) tibble())
  )
}

aa_to_1letter <- function(aa) {
  if (is.null(aa) || is.na(aa) || aa == "") return(NA_character_)
  aa <- as.character(aa)
  aa <- trimws(aa)
  
  if (nchar(aa) == 1 && toupper(aa) %in% AA_VALID) return(toupper(aa))
  
  if (nchar(aa) == 3) {
    if (aa %in% names(AA_3TO1)) return(unname(AA_3TO1[aa]))
    aa_title <- paste0(toupper(substr(aa, 1, 1)), tolower(substr(aa, 2, 3)))
    if (aa_title %in% names(AA_3TO1)) return(unname(AA_3TO1[aa_title]))
    if (toupper(aa) %in% names(AA_3TO1)) return(unname(AA_3TO1[toupper(aa)]))
  }
  NA_character_
}

parse_protein_change <- function(prot, title = NULL) {
  result <- list(position = NA_integer_, ref_aa = NA_character_, alt_aa = NA_character_)
  
  texts_to_try <- c(prot, title)
  
  for (txt in texts_to_try) {
    if (is.null(txt) || is.na(txt) || txt == "") next
    txt <- as.character(txt)
    
    m1 <- stringr::str_match(txt, "[\\(]?p\\.([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2})[\\)]?")
    if (!is.na(m1[1,1])) {
      result$position <- as.integer(m1[1,3])
      result$ref_aa <- aa_to_1letter(m1[1,2])
      result$alt_aa <- aa_to_1letter(m1[1,4])
      if (!is.na(result$ref_aa) && !is.na(result$alt_aa)) return(result)
    }
    
    m2 <- stringr::str_match(txt, "p\\.([A-Z])(\\d+)([A-Z])")
    if (!is.na(m2[1,1])) {
      result$position <- as.integer(m2[1,3])
      result$ref_aa <- m2[1,2]
      result$alt_aa <- m2[1,4]
      if (result$ref_aa %in% AA_VALID && result$alt_aa %in% AA_VALID) return(result)
    }
    
    m3 <- stringr::str_match(txt, "([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2})")
    if (!is.na(m3[1,1])) {
      result$position <- as.integer(m3[1,3])
      result$ref_aa <- aa_to_1letter(m3[1,2])
      result$alt_aa <- aa_to_1letter(m3[1,4])
      if (!is.na(result$ref_aa) && !is.na(result$alt_aa)) return(result)
    }
    
    m4 <- stringr::str_match(txt, "([A-Z])(\\d+)([A-Z])")
    if (!is.na(m4[1,1])) {
      pos_candidate <- as.integer(m4[1,3])
      ref_candidate <- m4[1,2]
      alt_candidate <- m4[1,4]
      if (ref_candidate %in% AA_VALID && alt_candidate %in% AA_VALID && pos_candidate > 0) {
        result$position <- pos_candidate
        result$ref_aa <- ref_candidate
        result$alt_aa <- alt_candidate
        return(result)
      }
    }
  }
  
  for (txt in texts_to_try) {
    if (is.null(txt) || is.na(txt) || txt == "") next
    nums <- stringr::str_extract_all(txt, "\\d+")[[1]]
    if (length(nums) > 0) {
      result$position <- as.integer(nums[1])
      break
    }
  }
  
  result
}

search_clinvar <- function(gene, max_vars = 100) {
  # Search for ALL variants (not just missense) to include pathogenic deletions, etc.
  url1 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=",
                 URLencode(paste0(gene, "[gene]")), 
                 "&retmode=json&retmax=", max_vars)
  r1 <- safe_get(url1)
  if (is.null(r1)) return(tibble())
  
  json1 <- jsonlite::fromJSON(httr::content(r1, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  ids <- json1$esearchresult$idlist
  
  if (length(ids) == 0) return(tibble())
  
  results <- list()
  
  for (id in ids) {
    Sys.sleep(0.15)  # Rate limiting
    url2 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", id, "&retmode=json")
    r2 <- safe_get(url2, timeout_secs = 10)
    if (is.null(r2)) next
    
    tryCatch({
      json2 <- jsonlite::fromJSON(httr::content(r2, "text", encoding = "UTF-8"), simplifyVector = FALSE)
      v <- json2$result[[id]]
      if (is.null(v)) next
      
      prot <- v$protein_change %||% ""
      title_str <- v$title %||% ""
      
      sig <- tryCatch({
        if (is.list(v$germline_classification)) v$germline_classification$description
        else as.character(v$germline_classification)
      }, error = function(e) NA)
      
      parsed <- parse_protein_change(prot, title_str)
      
      results[[length(results) + 1]] <- tibble(
        variation_id = as.character(id),
        title = as.character(title_str),
        protein_change = as.character(prot),
        significance = as.character(sig %||% "Not provided"),
        variant_type = as.character(v$obj_type %||% "Unknown"),
        position = parsed$position,
        ref_aa = parsed$ref_aa,
        alt_aa = parsed$alt_aa
      )
    }, error = function(e) NULL)
  }
  
  if (length(results) == 0) return(tibble())
  bind_rows(results)
}

# AlphaMissense with rate limiting and retry
get_alphamissense <- function(uniprot_acc, position, ref_aa, alt_aa, retry = 2) {
  if (is.null(position) || is.na(position)) {
    return(list(status = "no_data", message = "No position available"))
  }
  if (is.null(ref_aa) || is.na(ref_aa) || ref_aa == "") {
    return(list(status = "no_data", message = "Reference AA not available"))
  }
  if (is.null(alt_aa) || is.na(alt_aa) || alt_aa == "") {
    return(list(status = "no_data", message = "Alternate AA not available"))
  }
  
  ref_1 <- if (nchar(ref_aa) == 1) toupper(ref_aa) else aa_to_1letter(ref_aa)
  alt_1 <- if (nchar(alt_aa) == 1) toupper(alt_aa) else aa_to_1letter(alt_aa)
  
  if (is.na(ref_1) || !ref_1 %in% AA_VALID) {
    return(list(status = "no_data", message = paste("Invalid ref AA:", ref_aa)))
  }
  if (is.na(alt_1) || !alt_1 %in% AA_VALID) {
    return(list(status = "no_data", message = paste("Invalid alt AA:", alt_aa)))
  }
  
  hgvsp <- paste0(uniprot_acc, ":p.", ref_1, position, alt_1)
  url <- paste0("https://rest.ensembl.org/vep/human/hgvs/", URLencode(hgvsp), 
                "?content-type=application/json&AlphaMissense=1")
  
  for (attempt in 1:retry) {
    if (attempt > 1) Sys.sleep(2)  # Wait before retry
    
    r <- safe_get(url, timeout_secs = 30)
    if (is.null(r)) next
    
    tryCatch({
      txt <- httr::content(r, "text", encoding = "UTF-8")
      json <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
      
      if (!is.null(json$error)) {
        if (grepl("rate limit", tolower(json$error))) {
          Sys.sleep(3)
          next
        }
        return(list(status = "error", message = json$error))
      }
      
      if (length(json) > 0 && !is.null(json[[1]]$transcript_consequences)) {
        for (tc in json[[1]]$transcript_consequences) {
          if (!is.null(tc$alphamissense)) {
            am <- tc$alphamissense
            score <- am$am_pathogenicity
            cls <- am$am_class
            if (!is.null(score)) {
              return(list(
                status = "found",
                score = score,
                classification = cls,
                message = paste0(cls, " (", round(score, 3), ")")
              ))
            }
          }
        }
      }
      return(list(status = "not_found", message = "No AlphaMissense data for this variant"))
    }, error = function(e) NULL)
  }
  
  list(status = "error", message = "Ensembl API unavailable (rate limited). Try again in a few seconds.")
}

get_pdbe_mappings <- function(accession) {
  url <- glue::glue("https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/{accession}")
  r <- safe_get(url)
  if (is.null(r)) return(tibble())
  
  json <- tryCatch(
    jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(json)) return(tibble())
  
  rows <- list()
  for (acc in names(json)) {
    entry <- json[[acc]]
    if (!"PDB" %in% names(entry)) next
    
    for (pdb_id in names(entry$PDB)) {
      for (chain_data in entry$PDB[[pdb_id]]) {
        tryCatch({
          rows[[length(rows) + 1]] <- tibble(
            pdb = pdb_id,
            chain = chain_data$chain_id %||% "A",
            unp_start = as.integer(chain_data$unp_start %||% NA),
            unp_end = as.integer(chain_data$unp_end %||% NA),
            pdb_start = as.integer(chain_data$start$residue_number %||% chain_data$start %||% NA),
            identity = as.numeric(chain_data$identity %||% NA)
          )
        }, error = function(e) {})
      }
    }
  }
  
  if (length(rows) == 0) return(tibble())
  
  result <- bind_rows(rows)
  max_len <- max(result$unp_end, na.rm = TRUE)
  result %>%
    mutate(coverage = (unp_end - unp_start + 1) / max_len) %>%
    arrange(desc(coverage))
}

get_pdbe_mutations <- function(pdb_id) {
  url <- glue::glue("https://www.ebi.ac.uk/pdbe/api/pdb/entry/mutated_AA_or_nucleotide/{tolower(pdb_id)}")
  r <- safe_get(url, timeout_secs = 15)
  if (is.null(r)) return(NULL)
  
  tryCatch({
    txt <- httr::content(r, "text", encoding = "UTF-8")
    json <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
    json[[tolower(pdb_id)]]
  }, error = function(e) NULL)
}

get_pdb_summary <- function(pdb_id) {
  url <- glue::glue("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{tolower(pdb_id)}")
  r <- safe_get(url, timeout_secs = 15)
  if (is.null(r)) return(NULL)
  
  tryCatch({
    txt <- httr::content(r, "text", encoding = "UTF-8")
    json <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
    json[[tolower(pdb_id)]][[1]]
  }, error = function(e) NULL)
}

# ============================================================
# PRIORITY SCORE CALCULATION
# ============================================================
calculate_priority_score <- function(variants_df, domains_df, pdbs_df, plddt_data = NULL) {
  if (is.null(variants_df) || nrow(variants_df) == 0) return(variants_df)
  
  # Ensure variant_type column exists
  if (!"variant_type" %in% names(variants_df)) {
    variants_df$variant_type <- "Unknown"
  }
  
  variants_df <- variants_df %>%
    mutate(
      # Clinical severity score (0-40 points)
      clinical_score = case_when(
        grepl("^Pathogenic$", significance, ignore.case = TRUE) ~ 40,
        grepl("Pathogenic/Likely pathogenic", significance, ignore.case = TRUE) ~ 38,
        grepl("Likely pathogenic", significance, ignore.case = TRUE) ~ 32,
        grepl("risk factor", significance, ignore.case = TRUE) ~ 24,
        grepl("Uncertain", significance, ignore.case = TRUE) ~ 16,
        grepl("Conflicting", significance, ignore.case = TRUE) ~ 12,
        grepl("Likely benign", significance, ignore.case = TRUE) ~ 4,
        grepl("Benign/Likely benign", significance, ignore.case = TRUE) ~ 2,
        grepl("^Benign$", significance, ignore.case = TRUE) ~ 0,
        TRUE ~ 8
      ),
      
      # Domain membership score (0-25 points)
      in_domain = FALSE,
      domain_name = NA_character_
    )
  
  # Check domain membership
  if (!is.null(domains_df) && nrow(domains_df) > 0) {
    for (i in seq_len(nrow(variants_df))) {
      pos <- variants_df$position[i]
      if (is.na(pos)) next
      
      for (j in seq_len(nrow(domains_df))) {
        if (!is.na(domains_df$start[j]) && !is.na(domains_df$end[j])) {
          if (pos >= domains_df$start[j] && pos <= domains_df$end[j]) {
            variants_df$in_domain[i] <- TRUE
            variants_df$domain_name[i] <- domains_df$description[j]
            break
          }
        }
      }
    }
  }
  
  variants_df <- variants_df %>%
    mutate(domain_score = ifelse(in_domain, 25, 0))
  
  # PDB coverage score (0-15 points)
  variants_df$pdb_score <- 0
  if (!is.null(pdbs_df) && nrow(pdbs_df) > 0) {
    for (i in seq_len(nrow(variants_df))) {
      pos <- variants_df$position[i]
      if (is.na(pos)) next
      
      for (j in seq_len(nrow(pdbs_df))) {
        if (!is.na(pdbs_df$unp_start[j]) && !is.na(pdbs_df$unp_end[j])) {
          if (pos >= pdbs_df$unp_start[j] && pos <= pdbs_df$unp_end[j]) {
            variants_df$pdb_score[i] <- 15
            break
          }
        }
      }
    }
  }
  
  # Amino acid change severity (0-20 points)
  variants_df <- variants_df %>%
    mutate(
      aa_change_score = {
        scores <- numeric(n())
        for (i in seq_len(n())) {
          ref <- ref_aa[i]
          alt <- alt_aa[i]
          if (is.na(ref) || is.na(alt)) {
            scores[i] <- 10  # Unknown
          } else {
            # Calculate property differences
            hydro_diff <- abs(AA_PROPERTIES$hydrophobicity[ref] - AA_PROPERTIES$hydrophobicity[alt])
            size_diff <- abs(AA_PROPERTIES$size[ref] - AA_PROPERTIES$size[alt])
            charge_diff <- abs(AA_PROPERTIES$charge[ref] - AA_PROPERTIES$charge[alt])
            
            # Weighted sum (charge changes are most severe)
            scores[i] <- round((hydro_diff * 5 + size_diff * 5 + charge_diff * 10), 1)
          }
        }
        scores
      }
    )
  
  # Calculate total priority score
  variants_df <- variants_df %>%
    mutate(
      priority_score = clinical_score + domain_score + pdb_score + aa_change_score,
      priority_rank = case_when(
        priority_score >= 70 ~ "Critical",
        priority_score >= 50 ~ "High",
        priority_score >= 30 ~ "Medium",
        TRUE ~ "Low"
      )
    ) %>%
    arrange(desc(priority_score))
  
  variants_df
}

# ============================================================
# RAMACHANDRAN FUNCTIONS
# ============================================================
calculate_torsion_angles <- function(pdb_text, chain_filter = NULL) {
  tmpfile <- tempfile(fileext = ".pdb")
  writeLines(pdb_text, tmpfile)
  
  tryCatch({
    pdb <- bio3d::read.pdb(tmpfile, verbose = FALSE)
    
    chains <- unique(pdb$atom[, "chain"])
    if (!is.null(chain_filter) && chain_filter %in% chains) {
      chains <- chain_filter
    }
    
    torsion_all <- data.frame()
    
    for (ch in chains) {
      pdb_chain <- bio3d::trim.pdb(pdb, chain = ch)
      resno_aa <- unique(pdb_chain$atom[pdb_chain$atom[, "resid"] %in% ALL_AA, "resno"])
      if (length(resno_aa) == 0) next
      
      pdb_chain <- bio3d::trim.pdb(pdb_chain, resno = resno_aa)
      if (nrow(pdb_chain$atom) == 0) next
      
      tor <- bio3d::torsion.pdb(pdb_chain)
      tortab <- tor$tbl[, c("phi", "psi")]
      
      rn <- rownames(tortab)
      parsed <- strsplit(rn, ".", fixed = TRUE)
      
      df <- data.frame(
        phi = tortab[, "phi"],
        psi = tortab[, "psi"],
        resno = sapply(parsed, function(x) as.integer(x[1])),
        chain = sapply(parsed, function(x) x[2]),
        resid = sapply(parsed, function(x) x[3]),
        stringsAsFactors = FALSE
      )
      
      torsion_all <- rbind(torsion_all, df)
    }
    
    unlink(tmpfile)
    if (nrow(torsion_all) == 0) return(NULL)
    
    torsion_all <- torsion_all[!is.na(torsion_all$phi) & !is.na(torsion_all$psi), ]
    
    torsion_all$category <- "General"
    torsion_all$category[torsion_all$resid == "GLY"] <- "Glycine"
    torsion_all$category[torsion_all$resid == "PRO"] <- "Proline"
    
    for (i in 1:(nrow(torsion_all) - 1)) {
      if (torsion_all$resid[i + 1] == "PRO" && 
          torsion_all$chain[i] == torsion_all$chain[i + 1]) {
        torsion_all$category[i] <- "Pre-Proline"
      }
    }
    
    torsion_all$label <- paste0(torsion_all$resid, " ", torsion_all$resno, 
                                " (", torsion_all$chain, ")")
    torsion_all
    
  }, error = function(e) {
    unlink(tmpfile)
    NULL
  })
}

create_ramachandran_plotly <- function(torsion_data, highlight_resno = NULL, title = "Ramachandran Plot") {
  if (is.null(torsion_data) || nrow(torsion_data) == 0) return(NULL)
  
  cat_colors <- c("General" = "#3182bd", "Glycine" = "#e6550d", 
                  "Proline" = "#31a354", "Pre-Proline" = "#756bb1")
  
  p <- plot_ly(source = "rama_plot") %>%
    layout(
      title = list(text = title, x = 0.5),
      xaxis = list(title = "φ (°)", range = c(-180, 180), dtick = 60, zeroline = TRUE),
      yaxis = list(title = "ψ (°)", range = c(-180, 180), dtick = 60, zeroline = TRUE),
      shapes = list(
        list(type = "rect", x0 = -160, x1 = -20, y0 = -120, y1 = 50, fillcolor = "#74A9CF", opacity = 0.3, line = list(width = 0)),
        list(type = "rect", x0 = -180, x1 = -40, y0 = 80, y1 = 180, fillcolor = "#74A9CF", opacity = 0.3, line = list(width = 0)),
        list(type = "rect", x0 = -180, x1 = -40, y0 = -180, y1 = -120, fillcolor = "#74A9CF", opacity = 0.3, line = list(width = 0)),
        list(type = "rect", x0 = 20, x1 = 120, y0 = -60, y1 = 80, fillcolor = "#BDC9E1", opacity = 0.2, line = list(width = 0))
      ),
      annotations = list(
        list(x = -75, y = -40, text = "α", showarrow = FALSE, font = list(size = 20, color = "gray50")),
        list(x = -110, y = 140, text = "β", showarrow = FALSE, font = list(size = 20, color = "gray50"))
      ),
      plot_bgcolor = "#F8F9FA",
      margin = list(t = 50, b = 50, l = 60, r = 20)
    )
  
  for (cat in unique(torsion_data$category)) {
    df_cat <- torsion_data[torsion_data$category == cat, ]
    p <- p %>% add_trace(data = df_cat, x = ~phi, y = ~psi, type = "scatter", mode = "markers",
                         marker = list(color = cat_colors[cat], size = 8, opacity = 0.7), text = ~label,
                         hovertemplate = "<b>%{text}</b><br>φ: %{x:.1f}°<br>ψ: %{y:.1f}°<extra></extra>",
                         name = cat, key = ~resno)
  }
  
  if (!is.null(highlight_resno)) {
    hl <- torsion_data[torsion_data$resno == highlight_resno, ]
    if (nrow(hl) > 0) {
      p <- p %>% add_trace(data = hl, x = ~phi, y = ~psi, type = "scatter", mode = "markers",
                           marker = list(color = "red", size = 15, symbol = "circle-open", line = list(width = 3)),
                           text = ~label, hovertemplate = "<b>SELECTED: %{text}</b><extra></extra>",
                           name = "Selected", showlegend = FALSE)
    }
  }
  
  p %>% config(displayModeBar = TRUE, displaylogo = FALSE) %>% event_register("plotly_click")
}

# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .disclaimer-box { position: fixed; bottom: 10px; left: 10px; background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px 15px; max-width: 300px; font-size: 11px; color: #495057; z-index: 1000; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
      .plddt-legend { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px; padding: 8px; background-color: #f8f9fa; border-radius: 5px; font-size: 10px; }
      .plddt-item { display: flex; align-items: center; gap: 4px; }
      .plddt-color { width: 14px; height: 14px; border-radius: 3px; border: 1px solid #ccc; }
      .external-link { color: #007bff; cursor: pointer; text-decoration: none; }
      .external-link:hover { text-decoration: underline; }
      .mutation-warning { background-color: #fff3cd; border-left: 4px solid #ffc107; padding: 12px; margin: 10px 0; border-radius: 4px; font-size: 12px; }
      .assembly-info { background-color: #e3f2fd; border-left: 4px solid #2196f3; padding: 12px; margin: 10px 0; border-radius: 4px; font-size: 12px; }
      .view-options { background: #e9ecef; padding: 10px; border-radius: 5px; margin-top: 10px; }
      .tips-box { font-size: 11px; color: #666; background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px; }
      .variant-panel { background: #f0f7ff; padding: 12px; border-radius: 5px; margin-top: 10px; border-left: 3px solid #2196f3; }
      .links-panel { background: #f8f9fa; padding: 12px; border-radius: 5px; margin-top: 15px; }
      .am-benign { background-color: #c8e6c9; color: #2e7d32; padding: 3px 8px; border-radius: 3px; font-weight: bold; }
      .am-ambiguous { background-color: #fff9c4; color: #f9a825; padding: 3px 8px; border-radius: 3px; font-weight: bold; }
      .am-pathogenic { background-color: #ffcdd2; color: #c62828; padding: 3px 8px; border-radius: 3px; font-weight: bold; }
      .am-loading { color: #666; font-style: italic; }
      .variant-highlight-info { background: #e8f5e9; padding: 8px; border-radius: 4px; font-size: 11px; margin-top: 10px; border-left: 3px solid #4caf50; }
      .priority-critical { background-color: #d32f2f; color: white; padding: 2px 6px; border-radius: 3px; font-weight: bold; font-size: 11px; }
      .priority-high { background-color: #f57c00; color: white; padding: 2px 6px; border-radius: 3px; font-weight: bold; font-size: 11px; }
      .priority-medium { background-color: #fbc02d; color: black; padding: 2px 6px; border-radius: 3px; font-weight: bold; font-size: 11px; }
      .priority-low { background-color: #4caf50; color: white; padding: 2px 6px; border-radius: 3px; font-weight: bold; font-size: 11px; }
      .mutation-compare { background: #fff8e1; padding: 12px; border-radius: 5px; margin-top: 10px; border-left: 3px solid #ff9800; }
    "))
  ),
  
  titlePanel("ClinFold: Structural Impact Viewer for Genetic Variants"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. Protein Search"),
      textInput("gene", "Gene Symbol", value = "BRCA1", placeholder = "e.g., BRCA1, TP53"),
      tags$div(style = "font-size: 11px; color: #6c757d; margin-top: 5px;",
               tags$span(class = "external-link",
                         onclick = "window.open('https://www.genenames.org/about/guidelines/', '_blank')",
                         icon("external-link-alt"), " Gene nomenclature guidelines (HGNC)"
               )
      ),
      br(),
      actionButton("btn_search", "Search", class = "btn-primary", style = "width: 100%;"),
      hr(),
      
      uiOutput("protein_info_ui"),
      
      h4("2. Visualization"),
      checkboxInput("show_domains", "Show domains (colored)", value = TRUE),
      checkboxInput("show_ligands", "Show ligands (PDB only)", value = FALSE),
      checkboxInput("color_by_plddt", "Color by pLDDT (AlphaFold)", value = FALSE),
      
      conditionalPanel(
        condition = "input.color_by_plddt == true",
        tags$div(class = "plddt-legend",
                 tags$div(class = "plddt-item", tags$div(class = "plddt-color", style = "background-color: #0053D6;"), tags$span("Very high (>90)")),
                 tags$div(class = "plddt-item", tags$div(class = "plddt-color", style = "background-color: #65CBF3;"), tags$span("High (70-90)")),
                 tags$div(class = "plddt-item", tags$div(class = "plddt-color", style = "background-color: #FFDB13;"), tags$span("Low (50-70)")),
                 tags$div(class = "plddt-item", tags$div(class = "plddt-color", style = "background-color: #FF7D45;"), tags$span("Very low (<50)"))
        )
      ),
      
      tags$div(class = "view-options",
               h5("3D Animation"),
               radioButtons("animation_mode", NULL,
                            choices = c("Static" = "static", "Rocking" = "rocking", "Spinning" = "spinning"),
                            selected = "rocking", inline = TRUE)
      ),
      
      hr(),
      actionButton("btn_alphafold", "Load AlphaFold", class = "btn-success", style = "width: 100%;"),
      
      uiOutput("selected_variant_panel"),
      
      uiOutput("mutation_comparison_panel"),
      
      tags$div(class = "tips-box",
               tags$strong(icon("lightbulb"), " Tips:"),
               tags$ul(style = "padding-left: 15px; margin-top: 5px;",
                       tags$li("Search → AlphaFold → Select variant"),
                       tags$li("Click Ramachandran points → inspect residues"),
                       tags$li("Click a PDB → load experimental structure")
               )
      ),
      
      hr(),
      tags$div(style = "font-size: 10px; color: #666;",
               tags$strong("References:"),
               tags$ul(style = "padding-left: 15px; margin-top: 5px;",
                       tags$li(tags$span(class = "external-link", onclick = "window.open('https://github.com/BiKC/RamplotR', '_blank')", "RamplotR (BiKC)")),
                       tags$li(tags$span(class = "external-link", onclick = "window.open('https://naturegeorge.github.io/eigenblog/posts/introduce-pdb-profiling/', '_blank')", "pdb-profiling"))
               )
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel("Structural Overview", value = "tab_main",
                 br(),
                 uiOutput("structure_status"),
                 uiOutput("mutation_warning"),
                 uiOutput("assembly_info"),
                 uiOutput("variant_highlight_info"),
                 fluidRow(
                   column(6, h5("3D Structure"), NGLVieweR::NGLVieweROutput("ngl_view", height = "600px")),
                   column(6, h5("Ramachandran Plot"), plotlyOutput("rama_plot", height = "550px"), verbatimTextOutput("rama_click_info"))
                 )
        ),
        
        tabPanel("Clinical Variants", value = "tab_var",
                 br(), 
                 h4("ClinVar Variants with Priority Score"),
                 p("Variants are ranked by a ", tags$strong("Priority Score"), " combining clinical significance, domain location, PDB coverage, and amino acid change severity."),
                 DTOutput("variants_table"),
                 hr(),
                 uiOutput("priority_score_explanation"),
                 tags$div(class = "links-panel",
                          tags$a(href = "https://www.ncbi.nlm.nih.gov/clinvar/", target = "_blank", class = "btn btn-sm btn-outline-primary", icon("external-link-alt"), " ClinVar"),
                          tags$a(href = "https://alphamissense.hegelab.org/", target = "_blank", class = "btn btn-sm btn-outline-warning", icon("external-link-alt"), " AlphaMissense"),
                          tags$a(href = "https://gnomad.broadinstitute.org/", target = "_blank", class = "btn btn-sm btn-outline-info", icon("external-link-alt"), " gnomAD")
                 )
        ),
        
        tabPanel("PDB Structures", value = "tab_pdb",
                 br(), 
                 h4("Experimental Structures (PDBe/SIFTS)"),
                 p("Select a row to load the structure. ", 
                   tags$strong("Identity (%)"), " shows sequence match between UniProt and PDB."),
                 DTOutput("pdbs_table"),
                 hr(),
                 uiOutput("pdb_details_panel"),
                 tags$div(class = "links-panel",
                          tags$a(href = "https://www.rcsb.org/", target = "_blank", class = "btn btn-sm btn-outline-primary", icon("external-link-alt"), " RCSB PDB"),
                          tags$a(href = "https://www.ebi.ac.uk/pdbe/", target = "_blank", class = "btn btn-sm btn-outline-secondary", icon("external-link-alt"), " PDBe"),
                          tags$a(href = "https://alphafold.ebi.ac.uk/", target = "_blank", class = "btn btn-sm btn-outline-success", icon("external-link-alt"), " AlphaFold DB")
                 )
        ),
        
        tabPanel("UniProt Features", value = "tab_feat",
                 br(), 
                 h4("Protein Features from UniProt"),
                 DTOutput("features_table"),
                 tags$div(class = "links-panel",
                          uiOutput("uniprot_link"),
                          tags$a(href = "https://www.uniprot.org/", target = "_blank",
                                 class = "btn btn-sm btn-outline-primary", icon("external-link-alt"), " UniProt"),
                          tags$a(href = "https://www.ebi.ac.uk/interpro/", target = "_blank",
                                 class = "btn btn-sm btn-outline-secondary", icon("external-link-alt"), " InterPro")
                 )
        ),
        
        tabPanel("About", value = "tab_about",
                 br(),
                 h3("About ClinFold"),
                 tags$div(style = "max-width: 900px;",
                          p("ClinFold is a tool for visualizing the structural impact of genetic variants on protein structure. ",
                            "It integrates data from multiple authoritative databases to help researchers understand how mutations ",
                            "may affect protein function."),
                          
                          hr(),
                          h4(icon("calculator"), " Priority Score Calculation"),
                          p("The Priority Score (0-100) combines multiple evidence sources:"),
                          tags$table(class = "table table-bordered", style = "font-size: 12px;",
                                     tags$thead(tags$tr(tags$th("Component"), tags$th("Points"), tags$th("Description"))),
                                     tags$tbody(
                                       tags$tr(tags$td("Clinical Severity"), tags$td("0-40"), tags$td("Pathogenic=40, Likely pathogenic=32, Uncertain=16, Benign=0")),
                                       tags$tr(tags$td("Domain Location"), tags$td("0-25"), tags$td("25 points if variant is within a UniProt domain")),
                                       tags$tr(tags$td("PDB Coverage"), tags$td("0-15"), tags$td("15 points if variant position has experimental structure")),
                                       tags$tr(tags$td("AA Change Severity"), tags$td("0-20"), tags$td("Based on hydrophobicity, size, and charge differences"))
                                     )
                          ),
                          tags$ul(
                            tags$li(tags$span(class = "priority-critical", "Critical"), " - Score ≥ 70"),
                            tags$li(tags$span(class = "priority-high", "High"), " - Score 50-69"),
                            tags$li(tags$span(class = "priority-medium", "Medium"), " - Score 30-49"),
                            tags$li(tags$span(class = "priority-low", "Low"), " - Score < 30")
                          ),
                          
                          hr(),
                          h4(icon("dna"), " Mutation Comparison"),
                          p("When you select a variant, the sidebar shows a comparison of the wild-type and mutant amino acids, including:"),
                          tags$ul(
                            tags$li("Size difference (small → large or vice versa)"),
                            tags$li("Charge change (neutral → charged, etc.)"),
                            tags$li("Hydrophobicity change (polar → nonpolar)")
                          ),
                          p(tags$em("Note: AlphaFold cannot predict mutant structures in real-time. The 3D view shows the wild-type structure with the mutation site highlighted. ",
                                    "For mutant structure prediction, external tools like FoldX or ESMFold would be needed.")),
                          
                          hr(),
                          h4(icon("database"), " Data Sources"),
                          tags$ul(
                            tags$li(tags$strong("UniProt"), " - Curated protein sequences and functional annotations including domains, active sites, and binding regions."),
                            tags$li(tags$strong("ClinVar"), " - NCBI database of clinical significance of genetic variants with links to supporting evidence."),
                            tags$li(tags$strong("PDBe/SIFTS"), " - Experimental 3D structures from the Protein Data Bank with accurate sequence mappings via SIFTS."),
                            tags$li(tags$strong("AlphaFold"), " - AI-predicted protein structures from DeepMind, with per-residue confidence scores (pLDDT).")
                          ),
                          
                          hr(),
                          h4(icon("book"), " References"),
                          tags$ul(
                            tags$li("Jumper, J. et al. (2021). AlphaFold. Nature 596, 583–589."),
                            tags$li("Kyte, J. & Doolittle, R.F. (1982). A simple method for displaying the hydropathic character of a protein."),
                            tags$li("Zhu, Z. (2020). Introduction to pdb-profiling (blog post). EigenBlog. Available at https://naturegeorge.github.io/eigenblog/posts/introduce-pdb-profiling/"),
                            tags$li("BiKC (2025). RamplotR: R Shiny app for making Ramachandran plots. GitHub repository. Available at https://github.com/BiKC/RamplotR")
                          ),
                          
                          hr(),
                          h4(icon("user"), " Contact & Disclaimer"),
                          p("This tool is designed for ", tags$strong("research and educational purposes only"), ". ",
                            "Structural predictions and variant annotations are provided as reference information. ",
                            "Clinical decisions should always be based on validated diagnostic procedures and professional medical judgment."),
                          
                          br(),
                          
                          p(
                            tags$strong("Developer: "), "Fernando de la Puente Alonso de la Torre"
                          ),
                          p(
                            tags$strong("Email: "), 
                            tags$a(href = "mailto:puentealonsode@gmail.com", "puentealonsode@gmail.com")
                          ),
                          p(
                            tags$strong("Affiliation: "), "Cancer Epigenetics and Nanomedicine Laboratory"
                          )
                 )
        )
      )
    )
  ),
  
  tags$div(class = "disclaimer-box",
           tags$strong("Disclaimer"), tags$br(),
           "This tool is designed for research and educational purposes. ",
           "Structural predictions and variant annotations are provided as reference information only. ",
           "Clinical decisions should always be based on validated diagnostic procedures and professional medical judgment.",
           tags$br(), tags$br()
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    pdb_text = NULL, pdb_id = NULL, chain_id = NULL,
    unp_start = NULL, unp_end = NULL, pdb_start = NULL,
    identity = NULL, is_alphafold = FALSE,
    torsion_data = NULL, mutations = NULL, domains = NULL,
    pdb_summary = NULL, uniprot_acc = NULL,
    selected_rama_residue = NULL, 
    highlight_pdb_pos = NULL, highlight_chain = NULL,
    selected_variant_data = NULL, am_result = NULL, am_loading = FALSE,
    ngl_render_trigger = 0,
    variants_with_scores = NULL
  )
  
  # Search UniProt
  protein_data <- eventReactive(input$btn_search, {
    req(input$gene)
    rv$pdb_text <- NULL; rv$torsion_data <- NULL; rv$mutations <- NULL
    rv$selected_rama_residue <- NULL; rv$highlight_pdb_pos <- NULL
    rv$selected_variant_data <- NULL; rv$am_result <- NULL
    rv$highlight_chain <- NULL; rv$variants_with_scores <- NULL
    
    showNotification("Searching...", type = "message", duration = 2)
    result <- search_uniprot(input$gene)
    
    if (is.null(result)) {
      showNotification(paste("No protein found for", input$gene), type = "warning")
      return(NULL)
    }
    
    rv$uniprot_acc <- result$accession
    rv$domains <- if (!is.null(result$features) && nrow(result$features) > 0) {
      result$features %>% filter(type == "Domain")
    } else tibble()
    
    showNotification(paste("Found:", result$id), type = "message", duration = 2)
    result
  }, ignoreNULL = FALSE)
  
  # ClinVar variants - load automatically with search
  variants_data <- eventReactive(input$btn_search, {
    prot <- protein_data()
    if (is.null(prot)) return(tibble())
    
    showNotification("Loading ClinVar variants...", type = "message", duration = 2)
    vars <- withProgress(message = "Fetching variants...", { search_clinvar(prot$gene) })
    vars
  })
  
  # PDB structures - load automatically
  pdbs_data <- reactive({
    prot <- protein_data()
    req(prot)
    get_pdbe_mappings(prot$accession)
  })
  
  # Calculate priority scores when variants or PDBs change
  observe({
    vars <- variants_data()
    if (is.null(vars) || nrow(vars) == 0) {
      rv$variants_with_scores <- NULL
      return()
    }
    
    pdbs <- pdbs_data()
    domains <- rv$domains
    
    rv$variants_with_scores <- calculate_priority_score(vars, domains, pdbs)
  })
  
  output$protein_info_ui <- renderUI({
    prot <- protein_data()
    if (is.null(prot)) return(NULL)
    
    hgnc_url <- paste0("https://www.genenames.org/data/gene-symbol-report/#!/symbol/", toupper(prot$gene))
    
    tagList(
      h4("Selected Protein:"),
      tags$div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px;",
               tags$strong(prot$id), tags$br(),
               tags$small(prot$accession), tags$br(),
               tags$em(prot$name), tags$br(),
               tags$small(paste("Length:", prot$length, "aa")), tags$br(),
               tags$small("Gene: ",
                          tags$span(class = "external-link",
                                    onclick = paste0("window.open('", hgnc_url, "', '_blank')"),
                                    prot$gene, icon("external-link-alt", style = "font-size: 9px;")))
      )
    )
  })
  
  output$uniprot_link <- renderUI({
    prot <- protein_data()
    if (is.null(prot)) return(NULL)
    tags$a(href = paste0("https://www.uniprot.org/uniprotkb/", prot$accession), 
           target = "_blank", class = "btn btn-sm btn-outline-success", 
           style = "margin-right: 10px;",
           icon("external-link-alt"), paste(" View", prot$id, "on UniProt"))
  })
  
  # Variants table with priority score
  output$variants_table <- renderDT({
    df <- rv$variants_with_scores
    if (is.null(df) || nrow(df) == 0) {
      return(datatable(data.frame(Message = "No variants found. Click 'Search' first."), 
                       options = list(dom = 't'), rownames = FALSE))
    }
    
    df_display <- df %>%
      mutate(
        Priority = paste0(
          '<span class="priority-', tolower(priority_rank), '">', priority_rank, '</span>',
          '<br><small>', round(priority_score), ' pts</small>'
        ),
        ClinVar = paste0('<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/', variation_id, 
                         '/" target="_blank">', variation_id, '</a>'),
        Mutation = ifelse(protein_change == "" | is.na(protein_change), 
                          ifelse(!is.na(ref_aa) & !is.na(alt_aa), paste0(ref_aa, position, alt_aa), "—"),
                          protein_change),
        Type = ifelse("variant_type" %in% names(df), 
                      substr(variant_type, 1, 15), "—"),
        Position = ifelse(is.na(position), "—", position),
        Significance = ifelse(is.na(significance) | significance == "", "Not provided", significance),
        Domain = ifelse(in_domain, paste0("✓ ", substr(domain_name, 1, 15)), "—")
      ) %>%
      select(Priority, ClinVar, Mutation, Type, Position, Significance, Domain)
    
    datatable(df_display, selection = "single", rownames = FALSE, escape = FALSE,
              options = list(
                pageLength = 15, 
                scrollX = TRUE, 
                order = list(list(0, 'desc')),
                columnDefs = list(
                  list(width = '70px', targets = 0),
                  list(width = '60px', targets = 1),
                  list(width = '100px', targets = 2),
                  list(width = '80px', targets = 3),
                  list(width = '50px', targets = 4),
                  list(width = '110px', targets = 5),
                  list(width = '80px', targets = 6),
                  # Color the Significance column (index 5) using JS
                  list(
                    targets = 5,
                    createdCell = DT::JS("
                      function(td, cellData, rowData, row, col) {
                        var val = cellData.toLowerCase();
                        if (val.indexOf('pathogenic') !== -1 && val.indexOf('likely') === -1 && val.indexOf('benign') === -1) {
                          $(td).css({'background-color': '#ffcccc', 'font-weight': 'bold'});
                        } else if (val.indexOf('likely pathogenic') !== -1) {
                          $(td).css({'background-color': '#ffe0cc', 'font-weight': 'bold'});
                        } else if (val.indexOf('likely benign') !== -1) {
                          $(td).css({'background-color': '#e6ffcc'});
                        } else if (val.indexOf('benign') !== -1 && val.indexOf('likely') === -1 && val.indexOf('pathogenic') === -1) {
                          $(td).css({'background-color': '#ccffcc'});
                        } else if (val.indexOf('uncertain') !== -1) {
                          $(td).css({'background-color': '#f5f5f5'});
                        } else if (val.indexOf('conflicting') !== -1) {
                          $(td).css({'background-color': '#e6e6ff'});
                        } else if (val.indexOf('risk') !== -1) {
                          $(td).css({'background-color': '#fff3e0'});
                        }
                      }
                    ")
                  )
                )),
              filter = 'top')
  })
  
  output$priority_score_explanation <- renderUI({
    df <- rv$variants_with_scores
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    n_critical <- sum(df$priority_rank == "Critical")
    n_high <- sum(df$priority_rank == "High")
    
    tags$div(class = "alert alert-info",
             tags$strong(icon("chart-bar"), " Summary: "),
             n_critical, " Critical, ", n_high, " High priority variants out of ", nrow(df), " total.",
             tags$br(),
             tags$small("Priority Score = Clinical (40) + Domain (25) + PDB (15) + AA Change (20)")
    )
  })
  
  # Selected variant handling
  observeEvent(input$variants_table_rows_selected, {
    sel <- input$variants_table_rows_selected
    if (is.null(sel)) {
      rv$highlight_pdb_pos <- NULL
      rv$highlight_chain <- NULL
      rv$selected_variant_data <- NULL
      rv$am_result <- NULL
      return()
    }
    
    df <- rv$variants_with_scores
    if (is.null(df) || nrow(df) == 0) return()
    if (sel > nrow(df)) return()
    
    var_data <- df[sel, ]
    rv$selected_variant_data <- var_data
    rv$am_result <- NULL
    
    unp_pos <- var_data$position
    
    # Check if structure is loaded
    if (is.null(rv$pdb_text)) {
      showNotification("Loading AlphaFold structure...", type = "message", duration = 2)
      # Auto-load AlphaFold
      prot <- protein_data()
      if (!is.null(prot)) {
        tryCatch({
          api_url <- paste0("https://alphafold.ebi.ac.uk/api/prediction/", prot$accession)
          r1 <- safe_get(api_url)
          if (!is.null(r1)) {
            meta <- jsonlite::fromJSON(httr::content(r1, "text", encoding = "UTF-8"), simplifyVector = FALSE)
            pdb_url <- meta[[1]]$pdbUrl
            r2 <- safe_get(pdb_url, timeout_secs = 120)
            if (!is.null(r2)) {
              txt <- httr::content(r2, "text", encoding = "UTF-8")
              rv$pdb_text <- txt
              rv$pdb_id <- "AlphaFold"
              rv$chain_id <- "A"
              rv$unp_start <- 1
              rv$unp_end <- prot$length
              rv$pdb_start <- 1
              rv$identity <- 1.0
              rv$is_alphafold <- TRUE
              rv$mutations <- NULL
              rv$pdb_summary <- NULL
              rv$torsion_data <- calculate_torsion_angles(txt, "A")
            }
          }
        }, error = function(e) {
          showNotification("Could not load AlphaFold. Please load manually.", type = "warning")
        })
      }
    }
    
    if (is.na(unp_pos)) {
      showNotification("No position information for this variant", type = "warning")
      return()
    }
    
    pdb_pos <- NULL
    chain <- rv$chain_id
    
    if (!is.null(rv$pdb_text)) {
      if (!rv$is_alphafold && !is.null(rv$unp_start) && !is.null(rv$pdb_start)) {
        if (unp_pos >= rv$unp_start && unp_pos <= rv$unp_end) {
          pdb_pos <- rv$pdb_start + (unp_pos - rv$unp_start)
        } else {
          showNotification("Variant outside PDB coverage", type = "warning")
          return()
        }
      } else if (rv$is_alphafold) {
        pdb_pos <- unp_pos
        chain <- "A"
      }
      
      if (!is.null(pdb_pos)) {
        rv$highlight_pdb_pos <- pdb_pos
        rv$highlight_chain <- chain
        rv$ngl_render_trigger <- rv$ngl_render_trigger + 1
        
        updateTabsetPanel(session, "main_tabs", selected = "tab_main")
        showNotification(paste("Highlighting residue", pdb_pos), type = "message", duration = 2)
      }
    }
  })
  
  output$variant_highlight_info <- renderUI({
    if (is.null(rv$highlight_pdb_pos)) return(NULL)
    
    var_data <- rv$selected_variant_data
    ref <- if (!is.null(var_data)) var_data$ref_aa else "?"
    alt <- if (!is.null(var_data)) var_data$alt_aa else "?"
    
    tags$div(class = "variant-highlight-info",
             tags$strong(icon("crosshairs"), " Variant: "),
             ref, " → ", alt, " at position ", rv$highlight_pdb_pos, " (Chain ", rv$highlight_chain, ")",
             tags$br(),
             tags$small("Wild-type residue shown in element colors. Mutant structure cannot be predicted in real-time.")
    )
  })
  
  # Mutation comparison panel
  output$mutation_comparison_panel <- renderUI({
    var_data <- rv$selected_variant_data
    if (is.null(var_data)) return(NULL)
    
    ref <- var_data$ref_aa
    alt <- var_data$alt_aa
    
    if (is.na(ref) || is.na(alt)) return(NULL)
    
    # Get properties
    ref_hydro <- AA_PROPERTIES$hydrophobicity[ref]
    alt_hydro <- AA_PROPERTIES$hydrophobicity[alt]
    ref_size <- AA_PROPERTIES$size[ref]
    alt_size <- AA_PROPERTIES$size[alt]
    ref_charge <- AA_PROPERTIES$charge[ref]
    alt_charge <- AA_PROPERTIES$charge[alt]
    
    # Describe changes
    hydro_change <- if (abs(ref_hydro - alt_hydro) > 0.4) {
      if (ref_hydro > alt_hydro) "Hydrophobic → Polar" else "Polar → Hydrophobic"
    } else "Similar"
    
    size_change <- if (abs(ref_size - alt_size) > 0.4) {
      if (ref_size > alt_size) "Larger → Smaller" else "Smaller → Larger"
    } else "Similar"
    
    charge_change <- if (ref_charge != alt_charge) {
      paste0(ifelse(ref_charge > 0, "+", ifelse(ref_charge < 0, "-", "0")), " → ",
             ifelse(alt_charge > 0, "+", ifelse(alt_charge < 0, "-", "0")))
    } else "No change"
    
    tags$div(class = "mutation-compare",
             tags$strong(icon("exchange-alt"), " Mutation Effect"), tags$br(), tags$br(),
             tags$table(style = "font-size: 11px; width: 100%;",
                        tags$tr(tags$td(tags$strong("Property")), tags$td(tags$strong("Change"))),
                        tags$tr(tags$td("Hydrophobicity"), tags$td(hydro_change)),
                        tags$tr(tags$td("Size"), tags$td(size_change)),
                        tags$tr(tags$td("Charge"), tags$td(charge_change))
             ),
             tags$hr(),
             tags$small(tags$em("Note: 3D view shows wild-type. AlphaFold cannot predict mutant structures in real-time."))
    )
  })
  
  # Check AlphaMissense
  observeEvent(input$btn_check_am, {
    var_data <- rv$selected_variant_data
    if (is.null(var_data) || is.null(rv$uniprot_acc)) {
      showNotification("Please select a variant first", type = "warning")
      return()
    }
    
    rv$am_loading <- TRUE
    rv$am_result <- NULL
    
    result <- get_alphamissense(rv$uniprot_acc, var_data$position, var_data$ref_aa, var_data$alt_aa)
    
    rv$am_loading <- FALSE
    rv$am_result <- result
  })
  
  output$selected_variant_panel <- renderUI({
    var_data <- rv$selected_variant_data
    if (is.null(var_data)) return(NULL)
    
    prot_change <- var_data$protein_change
    pos <- var_data$position
    sig <- var_data$significance
    ref_aa <- var_data$ref_aa
    alt_aa <- var_data$alt_aa
    priority <- var_data$priority_rank
    score <- var_data$priority_score
    var_id <- var_data$variation_id
    var_type <- if ("variant_type" %in% names(var_data)) var_data$variant_type else "Unknown"
    
    has_aa <- !is.na(ref_aa) && !is.na(alt_aa)
    
    clinvar_url <- paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", var_id, "/")
    
    am_display <- NULL
    if (rv$am_loading) {
      am_display <- tags$span(class = "am-loading", icon("spinner", class = "fa-spin"), " Checking...")
    } else if (!is.null(rv$am_result)) {
      result <- rv$am_result
      if (result$status == "found") {
        css_class <- if (grepl("benign", tolower(result$classification))) "am-benign"
        else if (grepl("ambiguous", tolower(result$classification))) "am-ambiguous"
        else if (grepl("pathogenic", tolower(result$classification))) "am-pathogenic"
        else ""
        am_display <- tags$span(class = css_class, result$message)
      } else {
        am_display <- tags$span(style = "color: #666; font-size: 11px;", result$message)
      }
    }
    
    tags$div(class = "variant-panel",
             tags$strong(icon("dna"), " Selected Variant"), 
             tags$span(class = paste0("priority-", tolower(priority)), style = "margin-left: 10px;", priority),
             tags$br(), tags$br(),
             tags$span(style = "font-size: 14px;", 
                       ifelse(prot_change == "" || is.na(prot_change), 
                              ifelse(has_aa, paste0(ref_aa, pos, alt_aa), "No protein change"),
                              prot_change)), tags$br(),
             tags$small("Type: ", var_type), tags$br(),
             tags$small("Position: ", ifelse(is.na(pos), "N/A", pos), " | ", sig), tags$br(),
             tags$small("Score: ", round(score), " pts"), tags$br(), tags$br(),
             
             # ClinVar link button
             tags$a(href = clinvar_url, target = "_blank", 
                    class = "btn btn-sm btn-outline-primary", 
                    style = "margin-bottom: 8px; margin-right: 5px;",
                    icon("external-link-alt"), " View in ClinVar"),
             
             if (has_aa) {
               actionButton("btn_check_am", "AlphaMissense", class = "btn-sm btn-warning", 
                            icon = icon("search"), style = "margin-bottom: 8px;")
             } else {
               tags$span(style = "color: #999; font-size: 10px;", tags$br(), 
                         icon("info-circle"), " AlphaMissense: requires missense variant")
             },
             if (!is.null(am_display)) tags$div(style = "margin-top: 8px;", am_display)
    )
  })
  
  # Animation mode
  observeEvent(input$animation_mode, {
    if (is.null(rv$pdb_text)) return()
    
    if (input$animation_mode == "static") {
      NGLVieweR_proxy("ngl_view") %>% updateRock(FALSE) %>% updateSpin(FALSE)
    } else if (input$animation_mode == "rocking") {
      NGLVieweR_proxy("ngl_view") %>% updateSpin(FALSE) %>% updateRock(TRUE)
    } else if (input$animation_mode == "spinning") {
      NGLVieweR_proxy("ngl_view") %>% updateRock(FALSE) %>% updateSpin(TRUE)
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$show_ligands, {
    if (is.null(rv$pdb_text)) return()
    rv$ngl_render_trigger <- rv$ngl_render_trigger + 1
  }, ignoreInit = TRUE)
  
  output$pdbs_table <- renderDT({
    df <- pdbs_data()
    if (is.null(df) || nrow(df) == 0) {
      return(datatable(data.frame(Message = "No PDB structures"), options = list(dom = 't'), rownames = FALSE))
    }
    
    df %>%
      mutate(PDB = toupper(pdb), Chain = chain, 
             `UniProt Range` = paste(unp_start, "-", unp_end),
             `Identity (%)` = ifelse(!is.na(identity), round(identity * 100, 1), NA),
             `Coverage (%)` = round(coverage * 100, 1)) %>%
      select(PDB, Chain, `UniProt Range`, `Identity (%)`, `Coverage (%)`) %>%
      datatable(selection = "single", rownames = FALSE, 
                options = list(pageLength = 10, order = list(list(4, 'desc')))) %>%
      formatStyle('Identity (%)', backgroundColor = styleInterval(c(90, 99), c('#fff3cd', '#d4edda', '#d4edda')))
  })
  
  output$pdb_details_panel <- renderUI({
    sel <- input$pdbs_table_rows_selected
    if (is.null(sel)) return(tags$div(class = "alert alert-secondary", icon("hand-pointer"), " Select a PDB to load"))
    
    df <- pdbs_data()
    if (is.null(df) || nrow(df) == 0 || sel > nrow(df)) return(NULL)
    
    pdb_id <- df$pdb[sel]
    chain_id <- df$chain[sel]
    
    tags$div(class = "card", style = "padding: 15px; background: #e8f5e9;",
             tags$h5(icon("cube"), paste(" PDB", toupper(pdb_id), "- Chain", chain_id)),
             tags$a(href = paste0("https://www.rcsb.org/structure/", pdb_id), target = "_blank", 
                    class = "btn btn-sm btn-primary", icon("external-link-alt"), " RCSB")
    )
  })
  
  output$features_table <- renderDT({
    prot <- protein_data()
    if (is.null(prot) || is.null(prot$features) || nrow(prot$features) == 0) {
      return(datatable(data.frame(Message = "No features"), options = list(dom = 't'), rownames = FALSE))
    }
    
    feat <- prot$features
    if (!"description" %in% names(feat)) {
      feat$description <- NA_character_
    }
    
    feat %>% 
      filter(type != "Chain") %>%
      mutate(
        Type = type,
        Description = ifelse(is.na(description), "", description),
        Start = start,
        End = end
      ) %>%
      select(Type, Description, Start, End) %>%
      datatable(selection = "none", rownames = FALSE, options = list(pageLength = 15))
  })
  
  # Load PDB
  observeEvent(input$pdbs_table_rows_selected, {
    sel <- input$pdbs_table_rows_selected
    if (is.null(sel)) return()
    
    df <- pdbs_data()
    pdb_id <- df$pdb[sel]
    chain_id <- df$chain[sel]
    
    showNotification(paste("Loading", pdb_id, "..."), type = "message", duration = 2)
    
    tryCatch({
      url <- paste0("https://files.rcsb.org/download/", tolower(pdb_id), ".pdb")
      r <- safe_get(url, timeout_secs = 60)
      
      if (!is.null(r)) {
        txt <- httr::content(r, "text", encoding = "UTF-8")
        
        rv$pdb_text <- txt
        rv$pdb_id <- pdb_id
        rv$chain_id <- chain_id
        rv$unp_start <- df$unp_start[sel]
        rv$unp_end <- df$unp_end[sel]
        rv$pdb_start <- df$pdb_start[sel]
        rv$identity <- df$identity[sel]
        rv$is_alphafold <- FALSE
        rv$highlight_pdb_pos <- NULL
        rv$highlight_chain <- NULL
        
        rv$mutations <- get_pdbe_mutations(pdb_id)
        rv$pdb_summary <- get_pdb_summary(pdb_id)
        rv$torsion_data <- calculate_torsion_angles(txt, chain_id)
        
        updateTabsetPanel(session, "main_tabs", selected = "tab_main")
        showNotification(paste(pdb_id, "loaded"), type = "message", duration = 2)
      }
    }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
  })
  
  # Load AlphaFold
  observeEvent(input$btn_alphafold, {
    prot <- protein_data()
    req(prot)
    
    showNotification("Loading AlphaFold...", type = "message", duration = 2)
    
    tryCatch({
      api_url <- paste0("https://alphafold.ebi.ac.uk/api/prediction/", prot$accession)
      r1 <- safe_get(api_url)
      if (is.null(r1)) stop("AlphaFold not available")
      
      meta <- jsonlite::fromJSON(httr::content(r1, "text", encoding = "UTF-8"), simplifyVector = FALSE)
      pdb_url <- meta[[1]]$pdbUrl
      
      r2 <- safe_get(pdb_url, timeout_secs = 120)
      if (is.null(r2)) stop("Download failed")
      
      txt <- httr::content(r2, "text", encoding = "UTF-8")
      
      rv$pdb_text <- txt
      rv$pdb_id <- "AlphaFold"
      rv$chain_id <- "A"
      rv$unp_start <- 1
      rv$unp_end <- prot$length
      rv$pdb_start <- 1
      rv$identity <- 1.0
      rv$is_alphafold <- TRUE
      rv$mutations <- NULL
      rv$pdb_summary <- NULL
      rv$highlight_pdb_pos <- NULL
      rv$highlight_chain <- NULL
      
      rv$torsion_data <- calculate_torsion_angles(txt, "A")
      
      updateTabsetPanel(session, "main_tabs", selected = "tab_main")
      showNotification("AlphaFold loaded", type = "message", duration = 2)
    }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
  })
  
  output$structure_status <- renderUI({
    if (is.null(rv$pdb_text)) {
      return(tags$div(class = "alert alert-warning", icon("info-circle"), 
                      " Click 'Load AlphaFold' or select a PDB from the 'PDB Structures' tab"))
    }
    
    if (rv$is_alphafold) {
      return(tags$div(class = "alert alert-info", style = "margin-bottom: 10px;",
                      tags$strong(icon("brain"), " AlphaFold Structure"), tags$br(),
                      tags$small("Full protein coverage. Enable 'Color by pLDDT' to see confidence.")))
    }
    
    identity_str <- if (!is.null(rv$identity) && !is.na(rv$identity)) paste0(round(rv$identity * 100, 1), "%") else "N/A"
    
    tags$div(class = "alert alert-success", style = "margin-bottom: 10px;",
             tags$strong(icon("cube"), paste(" PDB:", toupper(rv$pdb_id), "(chain", rv$chain_id, ")")), tags$br(),
             tags$small(paste("UniProt:", rv$unp_start, "-", rv$unp_end, "| Identity:", identity_str)))
  })
  
  output$mutation_warning <- renderUI({
    if (is.null(rv$mutations) || length(rv$mutations) == 0) return(NULL)
    
    chain_muts <- list()
    for (mut in rv$mutations) {
      if (!is.null(mut$chain_id) && mut$chain_id == rv$chain_id) {
        chain_muts[[length(chain_muts) + 1]] <- mut
      }
    }
    if (length(chain_muts) == 0) return(NULL)
    
    mut_texts <- sapply(chain_muts, function(m) {
      paste0(m$mutation_details$from %||% "?", m$residue_number %||% "?", m$mutation_details$to %||% "?")
    })
    
    tags$div(class = "mutation-warning",
             tags$strong(icon("exclamation-triangle"), " ENGINEERED MUTATIONS"),
             tags$br(),
             tags$span(style = "font-size: 12px;", paste(mut_texts, collapse = ", ")))
  })
  
  output$assembly_info <- renderUI({
    if (is.null(rv$pdb_summary) || rv$is_alphafold) return(NULL)
    summ <- rv$pdb_summary
    
    tags$div(class = "assembly-info",
             tags$strong(icon("info-circle"), " "), tags$em(summ$title %||% ""))
  })
  
  # NGL Viewer
  output$ngl_view <- NGLVieweR::renderNGLVieweR({
    req(rv$pdb_text)
    rv$ngl_render_trigger
    
    tmpfile <- tempfile(fileext = ".pdb")
    writeLines(rv$pdb_text, tmpfile)
    
    ngl <- NGLVieweR(tmpfile) %>% stageParameters(backgroundColor = "white")
    
    if (rv$is_alphafold) {
      if (input$color_by_plddt) {
        ngl <- ngl %>% addRepresentation("cartoon", param = list(sele = "protein", colorScheme = "bfactor"))
      } else if (input$show_domains && !is.null(rv$domains) && nrow(rv$domains) > 0) {
        ngl <- ngl %>% addRepresentation("cartoon", param = list(sele = "protein", color = "lightgray", opacity = 0.5))
        for (i in seq_len(min(nrow(rv$domains), 10))) {
          dom <- rv$domains[i, ]
          if (!is.na(dom$start) && !is.na(dom$end)) {
            ngl <- ngl %>% addRepresentation("cartoon", param = list(
              sele = paste0(dom$start, "-", dom$end, ":A"), 
              color = DOMAIN_COLORS[(i %% length(DOMAIN_COLORS)) + 1], opacity = 1.0))
          }
        }
      } else {
        ngl <- ngl %>% addRepresentation("cartoon", param = list(sele = "protein", color = "spectrum"))
      }
    } else {
      if (input$show_domains && !is.null(rv$domains) && nrow(rv$domains) > 0) {
        ngl <- ngl %>%
          addRepresentation("cartoon", param = list(sele = paste0(":", rv$chain_id), color = "lightgray", opacity = 0.5)) %>%
          addRepresentation("cartoon", param = list(sele = paste0("not :", rv$chain_id), color = "lightgray", opacity = 0.15))
        
        for (i in seq_len(min(nrow(rv$domains), 10))) {
          dom <- rv$domains[i, ]
          if (!is.na(dom$start) && !is.na(dom$end) && dom$end >= rv$unp_start && dom$start <= rv$unp_end) {
            dom_pdb_start <- rv$pdb_start + max(0, dom$start - rv$unp_start)
            dom_pdb_end <- rv$pdb_start + min(dom$end - rv$unp_start, rv$unp_end - rv$unp_start)
            ngl <- ngl %>% addRepresentation("cartoon", param = list(
              sele = paste0(dom_pdb_start, "-", dom_pdb_end, ":", rv$chain_id),
              color = DOMAIN_COLORS[(i %% length(DOMAIN_COLORS)) + 1], opacity = 1.0))
          }
        }
      } else {
        ngl <- ngl %>%
          addRepresentation("cartoon", param = list(sele = paste0(":", rv$chain_id), color = "steelblue")) %>%
          addRepresentation("cartoon", param = list(sele = paste0("not :", rv$chain_id), color = "lightgray", opacity = 0.2))
      }
      
      if (!is.null(rv$mutations)) {
        for (mut in rv$mutations) {
          if (!is.null(mut$chain_id) && mut$chain_id == rv$chain_id && !is.null(mut$residue_number)) {
            ngl <- ngl %>% addRepresentation("ball+stick", param = list(
              sele = paste0(mut$residue_number, ":", rv$chain_id), color = "orange", radius = 0.3))
          }
        }
      }
      
      if (input$show_ligands) {
        ngl <- ngl %>% addRepresentation("ball+stick", param = list(sele = "ligand", colorScheme = "element"))
      }
    }
    
    # Variant highlight - THICK LICORICE
    if (!is.null(rv$highlight_pdb_pos) && !is.null(rv$highlight_chain)) {
      sele <- paste0(rv$highlight_pdb_pos, ":", rv$highlight_chain)
      ngl <- ngl %>% 
        addRepresentation("licorice", param = list(
          sele = sele, 
          colorScheme = "element",
          radius = 0.8,
          multipleBond = TRUE
        )) %>%
        addRepresentation("ball+stick", param = list(
          sele = sele,
          colorScheme = "element",
          aspectRatio = 1.5,
          bondScale = 0.4
        ))
      
      # Zoom to selection using zoomMove
      ngl <- ngl %>% zoomMove(sele, sele, 2000, 0)
    }
    
    if (input$animation_mode == "rocking") ngl <- ngl %>% setRock(TRUE)
    else if (input$animation_mode == "spinning") ngl <- ngl %>% setSpin(TRUE)
    
    ngl
  })
  
  # Ramachandran
  output$rama_plot <- renderPlotly({
    tor <- rv$torsion_data
    if (is.null(tor)) return(NULL)
    
    highlight <- rv$highlight_pdb_pos %||% rv$selected_rama_residue
    create_ramachandran_plotly(tor, highlight_resno = highlight, 
                               title = paste("Ramachandran -", toupper(rv$pdb_id %||% "")))
  })
  
  observeEvent(event_data("plotly_click", source = "rama_plot"), {
    click <- event_data("plotly_click", source = "rama_plot")
    if (!is.null(click$key)) rv$selected_rama_residue <- as.integer(click$key)
  })
  
  output$rama_click_info <- renderPrint({
    if (is.null(rv$selected_rama_residue)) { cat("Click a point for details"); return() }
    tor <- rv$torsion_data
    if (is.null(tor)) return()
    res <- tor[tor$resno == rv$selected_rama_residue, ]
    if (nrow(res) == 0) return()
    cat("Selected:", res$resid[1], res$resno[1], "(", res$chain[1], ")\n")
    cat("φ:", round(res$phi[1], 1), "° | ψ:", round(res$psi[1], 1), "°\n")
    cat("Category:", res$category[1])
  })
  
  outputOptions(output, "ngl_view", suspendWhenHidden = FALSE)
}

shinyApp(ui = ui, server = server)