# ============================================================
# ClinFold: Structural Impact Viewer for Genetic Variants
# Complete Shiny Application
# 
# References:
# - RamplotR (BiKC): https://github.com/BiKC/RamplotR
# - pdb-profiling: https://naturegeorge.github.io/eigenblog/posts/introduce-pdb-profiling/
# - Ramachandran et al. (1963)
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
AA_MAP <- c(
  ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLU="E", GLN="Q", GLY="G", HIS="H",
  ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W",
  TYR="Y", VAL="V"
)
ALL_AA <- names(AA_MAP)

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

search_clinvar <- function(gene, max_vars = 60) {
  url1 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=",
                 URLencode(paste0(gene, "[gene]")), "&retmode=json&retmax=", max_vars)
  r1 <- safe_get(url1)
  if (is.null(r1)) return(tibble())
  
  json1 <- jsonlite::fromJSON(httr::content(r1, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  ids <- json1$esearchresult$idlist
  if (length(ids) == 0) return(tibble())
  
  purrr::map_df(ids, function(id) {
    Sys.sleep(0.1)
    url2 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", id, "&retmode=json")
    r2 <- safe_get(url2, timeout_secs = 10)
    if (is.null(r2)) return(NULL)
    
    json2 <- jsonlite::fromJSON(httr::content(r2, "text", encoding = "UTF-8"), simplifyVector = FALSE)
    v <- json2$result[[id]]
    if (is.null(v)) return(NULL)
    
    prot <- v$protein_change %||% NA
    sig <- tryCatch({
      if (is.list(v$germline_classification)) v$germline_classification$description
      else as.character(v$germline_classification)
    }, error = function(e) NA)
    
    pos <- NA_integer_
    if (!is.na(prot) && prot != "") {
      nums <- stringr::str_extract_all(as.character(prot), "\\d+")[[1]]
      if (length(nums) > 0) pos <- as.integer(nums[1])
    }
    
    tibble(
      variation_id = id, 
      gene = gene,
      title = v$title %||% NA,
      protein_change = as.character(prot), 
      significance = sig, 
      position = pos
    )
  })
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
  
  cat_colors <- c(
    "General" = "#3182bd",
    "Glycine" = "#e6550d",
    "Proline" = "#31a354",
    "Pre-Proline" = "#756bb1"
  )
  
  p <- plot_ly(source = "rama_plot") %>%
    layout(
      title = list(text = title, x = 0.5),
      xaxis = list(title = "φ (°)", range = c(-180, 180), dtick = 60, zeroline = TRUE),
      yaxis = list(title = "ψ (°)", range = c(-180, 180), dtick = 60, zeroline = TRUE),
      shapes = list(
        list(type = "rect", x0 = -160, x1 = -20, y0 = -120, y1 = 50,
             fillcolor = "#74A9CF", opacity = 0.3, line = list(width = 0)),
        list(type = "rect", x0 = -180, x1 = -40, y0 = 80, y1 = 180,
             fillcolor = "#74A9CF", opacity = 0.3, line = list(width = 0)),
        list(type = "rect", x0 = -180, x1 = -40, y0 = -180, y1 = -120,
             fillcolor = "#74A9CF", opacity = 0.3, line = list(width = 0)),
        list(type = "rect", x0 = 20, x1 = 120, y0 = -60, y1 = 80,
             fillcolor = "#BDC9E1", opacity = 0.2, line = list(width = 0))
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
    p <- p %>% add_trace(
      data = df_cat, x = ~phi, y = ~psi,
      type = "scatter", mode = "markers",
      marker = list(color = cat_colors[cat], size = 8, opacity = 0.7),
      text = ~label,
      hovertemplate = paste0("<b>%{text}</b><br>φ: %{x:.1f}°<br>ψ: %{y:.1f}°<extra></extra>"),
      name = cat, key = ~resno
    )
  }
  
  if (!is.null(highlight_resno)) {
    hl <- torsion_data[torsion_data$resno == highlight_resno, ]
    if (nrow(hl) > 0) {
      p <- p %>% add_trace(
        data = hl, x = ~phi, y = ~psi,
        type = "scatter", mode = "markers",
        marker = list(color = "red", size = 15, symbol = "circle-open", line = list(width = 3)),
        text = ~label,
        hovertemplate = "<b>SELECTED: %{text}</b><br>φ: %{x:.1f}°<br>ψ: %{y:.1f}°<extra></extra>",
        name = "Selected", showlegend = FALSE
      )
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
      .disclaimer-box {
        position: fixed; bottom: 10px; left: 10px;
        background-color: #f8f9fa; border: 1px solid #dee2e6;
        border-radius: 5px; padding: 10px 15px; max-width: 300px;
        font-size: 11px; color: #495057; z-index: 1000;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      .plddt-legend {
        display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px;
        padding: 8px; background-color: #f8f9fa; border-radius: 5px; font-size: 10px;
      }
      .plddt-item { display: flex; align-items: center; gap: 4px; }
      .plddt-color { width: 14px; height: 14px; border-radius: 3px; border: 1px solid #ccc; }
      .external-link { color: #007bff; cursor: pointer; text-decoration: none; }
      .external-link:hover { text-decoration: underline; }
      .mutation-warning {
        background-color: #fff3cd; border-left: 4px solid #ffc107;
        padding: 12px; margin: 10px 0; border-radius: 4px; font-size: 12px;
      }
      .assembly-info {
        background-color: #e3f2fd; border-left: 4px solid #2196f3;
        padding: 12px; margin: 10px 0; border-radius: 4px; font-size: 12px;
      }
      .view-options { background: #e9ecef; padding: 10px; border-radius: 5px; margin-top: 10px; }
      .tips-box { font-size: 11px; color: #666; background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px; }
      .variant-highlight-info { background: #e8f5e9; padding: 8px; border-radius: 4px; font-size: 11px; margin-top: 10px; border-left: 3px solid #4caf50; }
      .links-panel { background: #f8f9fa; padding: 12px; border-radius: 5px; margin-top: 15px; }
      .links-panel a { margin-right: 15px; }
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
      
      h4("2. Visualization Options"),
      checkboxInput("show_domains", "Show domains (colored)", value = TRUE),
      checkboxInput("show_ligands", "Show ligands", value = FALSE),
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
               h5("3D View Options"),
               fluidRow(
                 column(6, checkboxInput("spinning", "Spinning", value = FALSE)),
                 column(6, checkboxInput("rocking", "Rocking", value = TRUE))
               )
      ),
      
      hr(),
      actionButton("btn_alphafold", "Load AlphaFold", class = "btn-success", style = "width: 100%;"),
      
      uiOutput("selected_variant_info"),
      
      tags$div(class = "tips-box",
               tags$strong(icon("info-circle"), " Tips:"),
               tags$ul(style = "padding-left: 15px; margin-top: 5px;",
                       tags$li("Click a variant → shown as sticks"),
                       tags$li("Click a PDB → load experimental structure"),
                       tags$li("Click Ramachandran points → inspect residues"),
                       tags$li(tags$span(style = "color: orange;", "ORANGE"), " = engineered mutations"),
                       tags$li("Domains colored, non-domain in gray"),
                       tags$li("Check mutation warnings for wild-type status")
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
        
        tabPanel("Structure & Ramachandran", value = "tab_main",
                 br(),
                 uiOutput("structure_status"),
                 uiOutput("mutation_warning"),
                 uiOutput("assembly_info"),
                 fluidRow(
                   column(6,
                          h5("3D Structure"),
                          NGLVieweR::NGLVieweROutput("ngl_view", height = "600px")
                   ),
                   column(6,
                          h5("Ramachandran Plot"),
                          plotlyOutput("rama_plot", height = "550px"),
                          verbatimTextOutput("rama_click_info")
                   )
                 )
        ),
        
        tabPanel("Clinical Variants", value = "tab_var",
                 br(), 
                 h4("ClinVar Variants"),
                 p("Select a row to highlight the variant position in the 3D structure (shown as sticks with element colors)."),
                 DTOutput("variants_table"),
                 tags$div(class = "links-panel",
                          tags$strong(icon("link"), " External Resources:"),
                          tags$br(), tags$br(),
                          tags$a(href = "https://www.ncbi.nlm.nih.gov/clinvar/", target = "_blank", 
                                 class = "btn btn-sm btn-outline-primary", icon("external-link-alt"), " ClinVar Database"),
                          tags$a(href = "https://www.ensembl.org/", target = "_blank",
                                 class = "btn btn-sm btn-outline-secondary", icon("external-link-alt"), " Ensembl"),
                          tags$a(href = "https://gnomad.broadinstitute.org/", target = "_blank",
                                 class = "btn btn-sm btn-outline-info", icon("external-link-alt"), " gnomAD")
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
                          tags$strong(icon("link"), " External Resources:"),
                          tags$br(), tags$br(),
                          tags$a(href = "https://www.rcsb.org/", target = "_blank",
                                 class = "btn btn-sm btn-outline-primary", icon("external-link-alt"), " RCSB PDB"),
                          tags$a(href = "https://www.ebi.ac.uk/pdbe/", target = "_blank",
                                 class = "btn btn-sm btn-outline-secondary", icon("external-link-alt"), " PDBe"),
                          tags$a(href = "https://alphafold.ebi.ac.uk/", target = "_blank",
                                 class = "btn btn-sm btn-outline-success", icon("external-link-alt"), " AlphaFold DB")
                 )
        ),
        
        tabPanel("UniProt Features", value = "tab_feat",
                 br(), 
                 h4("Protein Features from UniProt"),
                 p("Domains, binding sites, active sites, and other annotated regions from UniProt."),
                 DTOutput("features_table"),
                 tags$div(class = "links-panel",
                          tags$strong(icon("link"), " External Resources:"),
                          tags$br(), tags$br(),
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
                          
                          h4(icon("database"), " Data Sources"),
                          tags$ul(
                            tags$li(tags$strong("UniProt"), " - Curated protein sequences and functional annotations including domains, active sites, and binding regions."),
                            tags$li(tags$strong("ClinVar"), " - NCBI database of clinical significance of genetic variants with links to supporting evidence."),
                            tags$li(tags$strong("PDBe/SIFTS"), " - Experimental 3D structures from the Protein Data Bank with accurate sequence mappings via SIFTS."),
                            tags$li(tags$strong("AlphaFold"), " - AI-predicted protein structures from DeepMind, with per-residue confidence scores (pLDDT).")
                          ),
                          
                          hr(),
                          h4(icon("exclamation-triangle"), " Understanding PDB-UniProt Mappings"),
                          p("When mapping clinical variants to experimental structures, several important considerations apply:"),
                          
                          tags$div(style = "background: #fff3e0; padding: 15px; border-radius: 5px; margin: 10px 0;",
                                   h5("1. Partial Coverage"),
                                   p("PDB structures often cover only a fragment of the full protein. A variant may fall outside the crystallized region. ",
                                     "Always check the 'Coverage (%)' column to understand what portion of the protein is resolved.")
                          ),
                          
                          tags$div(style = "background: #ffebee; padding: 15px; border-radius: 5px; margin: 10px 0;",
                                   h5("2. Engineered Mutations"),
                                   p("Crystallographers frequently introduce mutations to improve protein stability, solubility, or crystallization. ",
                                     "These are shown in ", tags$span(style = "color: orange; font-weight: bold;", "ORANGE"), " in the 3D viewer."),
                                   tags$ul(
                                     tags$li("If a PDB has engineered mutations, it may NOT represent the wild-type structure"),
                                     tags$li("The 'Identity (%)' column indicates sequence match - values below 100% suggest modifications"),
                                     tags$li("A yellow warning panel appears when engineered mutations are detected in your selected chain")
                                   )
                          ),
                          
                          tags$div(style = "background: #e3f2fd; padding: 15px; border-radius: 5px; margin: 10px 0;",
                                   h5("3. Biological Assembly vs Asymmetric Unit"),
                                   p("The ", tags$strong("asymmetric unit"), " is the minimal repeating unit in the crystal. ",
                                     "The ", tags$strong("biological assembly"), " represents the functional oligomeric state of the protein."),
                                   tags$ul(
                                     tags$li("A protein may function as a dimer, but the asymmetric unit contains only one chain"),
                                     tags$li("Variants at protein-protein interfaces may disrupt oligomerization"),
                                     tags$li("Check the assembly information to understand the functional context")
                                   )
                          ),
                          
                          tags$div(style = "background: #f3e5f5; padding: 15px; border-radius: 5px; margin: 10px 0;",
                                   h5("4. Sequence Gaps"),
                                   p("Flexible loops or disordered regions may not be resolved in the electron density:"),
                                   tags$ul(
                                     tags$li("Missing residues appear as gaps in the structure"),
                                     tags$li("These regions may still be functionally important"),
                                     tags$li("AlphaFold can model these regions (check pLDDT for confidence)")
                                   )
                          ),
                          
                          hr(),
                          h4(icon("chart-area"), " Ramachandran Plot"),
                          p("The Ramachandran plot shows backbone torsion angles (φ, ψ) for each residue, ",
                            "first proposed by G.N. Ramachandran in 1963. It is used for:"),
                          tags$ul(
                            tags$li("Structure validation - most residues should fall in allowed regions"),
                            tags$li("Identifying unusual conformations that may indicate errors or functional importance"),
                            tags$li("Understanding how mutations might affect local backbone geometry")
                          ),
                          p("Different amino acids have characteristic distributions:"),
                          tags$ul(
                            tags$li(tags$span(style = "color: #e6550d;", "●"), " ", tags$strong("Glycine"), " - Most flexible, can adopt many conformations"),
                            tags$li(tags$span(style = "color: #31a354;", "●"), " ", tags$strong("Proline"), " - Most restricted due to cyclic side chain"),
                            tags$li(tags$span(style = "color: #756bb1;", "●"), " ", tags$strong("Pre-Proline"), " - Residues before proline have restricted φ angles")
                          ),
                          
                          hr(),
                          h4(icon("palette"), " pLDDT Confidence Score (AlphaFold)"),
                          p("The predicted Local Distance Difference Test (pLDDT) indicates per-residue prediction confidence:"),
                          tags$ul(
                            tags$li(tags$span(style = "color: #0053D6;", "■"), " ", tags$strong("Very high (>90)"), " - High confidence, backbone likely accurate"),
                            tags$li(tags$span(style = "color: #65CBF3;", "■"), " ", tags$strong("High (70-90)"), " - Good prediction, some uncertainty in side chains"),
                            tags$li(tags$span(style = "color: #FFDB13;", "■"), " ", tags$strong("Low (50-70)"), " - Caution, may be flexible or uncertain"),
                            tags$li(tags$span(style = "color: #FF7D45;", "■"), " ", tags$strong("Very low (<50)"), " - Likely disordered or unstructured")
                          ),
                          
                          hr(),
                          h4(icon("book"), " References"),
                          tags$ul(
                            tags$li("Ramachandran, G.N., Ramakrishnan, C. & Sasisekharan, V. (1963). ", 
                                    tags$em("Stereochemistry of polypeptide chain configurations."), " J. Mol. Biol. 7:95-99."),
                            tags$li(tags$a(href = "https://github.com/BiKC/RamplotR", target = "_blank", "RamplotR (BiKC)"), 
                                    " - Shiny app for Ramachandran plot visualization"),
                            tags$li(tags$a(href = "https://naturegeorge.github.io/eigenblog/posts/introduce-pdb-profiling/", target = "_blank", "pdb-profiling"), 
                                    " - Methodology for PDB structure profiling"),
                            tags$li(tags$a(href = "https://www.ebi.ac.uk/pdbe/docs/sifts/", target = "_blank", "SIFTS"), 
                                    " - Structure Integration with Function, Taxonomy and Sequences")
                          ),
                          
                          hr(),
                          h4(icon("user"), " Contact & Disclaimer"),
                          p("This tool is designed for ", tags$strong("research and educational purposes only"), ". ",
                            "Structural predictions and variant annotations are provided as reference information. ",
                            "Clinical decisions should always be based on validated diagnostic procedures and professional medical judgment.")
                 )
        )
      )
    )
  ),
  
  tags$div(class = "disclaimer-box",
           tags$strong("Academic Use Only"), tags$br(),
           "This tool is designed for research and educational purposes. ",
           "Structural predictions and variant annotations are provided as reference information only. ",
           "Clinical decisions should always be based on validated diagnostic procedures and professional medical judgment.",
           tags$br(), tags$br(),
           tags$small(style = "color: #6c757d;", "Data sources: UniProt, ClinVar, PDBe, AlphaFold")
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
    pdb_summary = NULL,
    selected_rama_residue = NULL, selected_variant_pos = NULL
  )
  
  # Search UniProt
  protein_data <- eventReactive(input$btn_search, {
    req(input$gene)
    rv$pdb_text <- NULL; rv$torsion_data <- NULL; rv$mutations <- NULL
    rv$selected_rama_residue <- NULL; rv$selected_variant_pos <- NULL
    rv$pdb_summary <- NULL
    
    showNotification("Searching UniProt...", type = "message", duration = 2)
    result <- search_uniprot(input$gene)
    
    if (is.null(result)) {
      showNotification(paste("No protein found for", input$gene), type = "warning")
      return(NULL)
    }
    
    if (!is.null(result$features) && nrow(result$features) > 0) {
      rv$domains <- result$features %>% filter(type == "Domain")
    } else {
      rv$domains <- tibble()
    }
    
    showNotification(paste("Found:", result$id), type = "message", duration = 2)
    result
  }, ignoreNULL = FALSE)
  
  # Protein info
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
    
    url <- paste0("https://www.uniprot.org/uniprotkb/", prot$accession)
    tags$a(href = url, target = "_blank", class = "btn btn-sm btn-outline-success", 
           style = "margin-right: 10px;", icon("external-link-alt"), paste(" View", prot$id, "on UniProt"))
  })
  
  # ClinVar variants
  variants_data <- eventReactive(input$btn_search, {
    prot <- protein_data()
    if (is.null(prot)) return(tibble())
    
    showNotification("Loading ClinVar...", type = "message", duration = 2)
    withProgress(message = "Fetching variants...", { df <- search_clinvar(prot$gene) })
    
    if (nrow(df) > 0) showNotification(paste("Found", nrow(df), "variants"), type = "message", duration = 2)
    df
  })
  
  # Variants table with adjusted column widths
  output$variants_table <- renderDT({
    df <- variants_data()
    if (nrow(df) == 0) {
      return(datatable(data.frame(Message = "No variants found"), options = list(dom = 't'), rownames = FALSE))
    }
    
    df_display <- df %>%
      mutate(
        ClinVar = paste0('<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/', variation_id, 
                         '/" target="_blank" title="View on ClinVar">', variation_id, 
                         ' <i class="fa fa-external-link-alt"></i></a>'),
        Mutation = protein_change,
        Position = position,
        Significance = significance,
        Title = title
      ) %>%
      select(ClinVar, Mutation, Position, Significance, Title) %>%
      arrange(Position)
    
    datatable(df_display, selection = "single", rownames = FALSE, escape = FALSE,
              options = list(
                pageLength = 12, 
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = list(
                  list(width = '70px', targets = 0),   # ClinVar
                  list(width = '150px', targets = 1),  # Mutation (WIDER)
                  list(width = '60px', targets = 2),   # Position
                  list(width = '100px', targets = 3),  # Significance (narrower)
                  list(width = '150px', targets = 4)   # Title (narrower)
                )
              ),
              filter = 'top') %>%
      formatStyle('Significance',
                  backgroundColor = styleEqual(
                    c('Pathogenic', 'Likely pathogenic', 'Benign', 'Likely benign', 'Uncertain significance'),
                    c('#ffcccc', '#ffe6cc', '#ccffcc', '#e6ffcc', '#f5f5f5')))
  })
  
  # Track selected variant - show with ELEMENT COLORS (like 3Dmol)
  observeEvent(input$variants_table_rows_selected, {
    sel <- input$variants_table_rows_selected
    if (is.null(sel)) { rv$selected_variant_pos <- NULL; return() }
    
    df <- variants_data()
    if (nrow(df) == 0) return()
    
    df_sorted <- df %>% arrange(position)
    if (sel > nrow(df_sorted)) return()
    
    unp_pos <- df_sorted$position[sel]
    
    if (!is.na(unp_pos) && !is.null(rv$pdb_text)) {
      if (!rv$is_alphafold && !is.null(rv$unp_start) && !is.null(rv$pdb_start)) {
        if (unp_pos >= rv$unp_start && unp_pos <= rv$unp_end) {
          pdb_pos <- rv$pdb_start + (unp_pos - rv$unp_start)
          rv$selected_variant_pos <- pdb_pos
          
          # Show as STICKS with ELEMENT colors (CPK-like coloring)
          NGLVieweR_proxy("ngl_view") %>%
            removeSelection("variant_highlight") %>%
            addSelection("licorice", 
                         param = list(
                           name = "variant_highlight",
                           sele = paste0(pdb_pos, ":", rv$chain_id),
                           colorScheme = "element",
                           radius = 0.3
                         ))
          
          showNotification(paste("Variant at PDB position", pdb_pos), type = "message", duration = 2)
        } else {
          showNotification("Variant outside PDB coverage", type = "warning")
          rv$selected_variant_pos <- NULL
        }
      } else if (rv$is_alphafold) {
        rv$selected_variant_pos <- unp_pos
        
        NGLVieweR_proxy("ngl_view") %>%
          removeSelection("variant_highlight") %>%
          addSelection("licorice",
                       param = list(
                         name = "variant_highlight", 
                         sele = paste0(unp_pos, ":A"),
                         colorScheme = "element",
                         radius = 0.3
                       ))
        
        showNotification(paste("Variant at position", unp_pos), type = "message", duration = 2)
      }
    }
  })
  
  output$selected_variant_info <- renderUI({
    if (is.null(rv$selected_variant_pos)) return(NULL)
    
    tags$div(class = "variant-highlight-info",
             tags$strong(icon("map-marker-alt"), " Variant Highlighted:"),
             tags$br(), paste("PDB residue:", rv$selected_variant_pos),
             tags$br(), tags$small("Shown as sticks with element colors")
    )
  })
  
  # Ligands toggle
  observeEvent(input$show_ligands, {
    if (is.null(rv$pdb_text)) return()
    
    if (input$show_ligands) {
      NGLVieweR_proxy("ngl_view") %>%
        addSelection("ball+stick", 
                     param = list(name = "ligands", sele = "ligand", colorScheme = "element"))
    } else {
      NGLVieweR_proxy("ngl_view") %>%
        removeSelection("ligands")
    }
  }, ignoreInit = TRUE)
  
  # PDB structures
  pdbs_data <- reactive({
    prot <- protein_data()
    req(prot)
    get_pdbe_mappings(prot$accession)
  })
  
  output$pdbs_table <- renderDT({
    df <- pdbs_data()
    if (nrow(df) == 0) {
      return(datatable(data.frame(Message = "No PDB structures found"), options = list(dom = 't'), rownames = FALSE))
    }
    
    df %>%
      mutate(
        PDB = toupper(pdb), Chain = chain,
        `UniProt Range` = paste(unp_start, "-", unp_end),
        `Identity (%)` = ifelse(!is.na(identity), round(identity * 100, 1), NA),
        `Coverage (%)` = round(coverage * 100, 1)
      ) %>%
      select(PDB, Chain, `UniProt Range`, `Identity (%)`, `Coverage (%)`) %>%
      datatable(selection = "single", rownames = FALSE, 
                options = list(pageLength = 10, order = list(list(4, 'desc')))) %>%
      formatStyle('Identity (%)',
                  backgroundColor = styleInterval(c(90, 99), c('#fff3cd', '#d4edda', '#d4edda')))
  })
  
  output$pdb_details_panel <- renderUI({
    sel <- input$pdbs_table_rows_selected
    if (is.null(sel)) {
      return(tags$div(class = "alert alert-secondary",
                      icon("hand-pointer"), " Select a PDB to see details and load it."))
    }
    
    df <- pdbs_data()
    if (nrow(df) == 0 || sel > nrow(df)) return(NULL)
    
    pdb_id <- df$pdb[sel]
    chain_id <- df$chain[sel]
    
    rcsb_url <- paste0("https://www.rcsb.org/structure/", pdb_id)
    pdbe_url <- paste0("https://www.ebi.ac.uk/pdbe/entry/pdb/", pdb_id)
    pdbe_3d <- paste0("https://www.ebi.ac.uk/pdbe/entry/view3D/", pdb_id, "/?viewer=ngl")
    
    tags$div(class = "card", style = "padding: 15px; background: #e8f5e9;",
             tags$h5(icon("cube"), paste(" PDB", toupper(pdb_id), "- Chain", chain_id)),
             tags$hr(),
             fluidRow(
               column(6,
                      tags$p(
                        tags$strong("UniProt Range: "), df$unp_start[sel], " - ", df$unp_end[sel], tags$br(),
                        tags$strong("Sequence Identity: "), 
                        if (!is.na(df$identity[sel])) paste0(round(df$identity[sel] * 100, 1), "%") else "N/A", tags$br(),
                        tags$strong("Coverage: "), round(df$coverage[sel] * 100, 1), "%"
                      )
               ),
               column(6,
                      tags$strong("Quick Links:"), tags$br(),
                      tags$a(href = rcsb_url, target = "_blank", class = "btn btn-sm btn-primary", style = "margin: 2px;", icon("external-link-alt"), " RCSB"),
                      tags$a(href = pdbe_url, target = "_blank", class = "btn btn-sm btn-secondary", style = "margin: 2px;", icon("external-link-alt"), " PDBe"),
                      tags$a(href = pdbe_3d, target = "_blank", class = "btn btn-sm btn-info", style = "margin: 2px;", icon("eye"), " 3D View")
               )
             ),
             tags$hr(),
             tags$small(style = "color: #666;",
                        icon("info-circle"), " Identity < 100% may indicate engineered mutations.")
    )
  })
  
  # Features table - NO ROW SELECTION
  output$features_table <- renderDT({
    prot <- protein_data()
    if (is.null(prot) || is.null(prot$features) || nrow(prot$features) == 0) {
      return(datatable(data.frame(Message = "No features found"), options = list(dom = 't'), rownames = FALSE))
    }
    
    prot$features %>%
      filter(type != "Chain") %>%
      select(type, description, start, end) %>%
      rename(Type = type, Description = description, Start = start, End = end) %>%
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
        rv$selected_variant_pos <- NULL
        
        rv$mutations <- get_pdbe_mutations(pdb_id)
        rv$pdb_summary <- get_pdb_summary(pdb_id)
        rv$torsion_data <- calculate_torsion_angles(txt, chain_id)
        
        showNotification(paste(pdb_id, "loaded"), type = "message", duration = 2)
        updateTabsetPanel(session, "main_tabs", selected = "tab_main")
      }
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
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
      if (is.null(pdb_url)) stop("No PDB URL")
      
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
      rv$selected_variant_pos <- NULL
      
      rv$torsion_data <- calculate_torsion_angles(txt, "A")
      
      showNotification("AlphaFold loaded", type = "message", duration = 2)
      updateTabsetPanel(session, "main_tabs", selected = "tab_main")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Structure status
  output$structure_status <- renderUI({
    if (is.null(rv$pdb_text)) {
      return(tags$div(class = "alert alert-warning",
                      icon("info-circle"), " Click 'Load AlphaFold' or select a PDB from the 'PDB Structures' tab"))
    }
    
    if (rv$is_alphafold) {
      return(tags$div(class = "alert alert-info", style = "margin-bottom: 10px;",
                      tags$strong(icon("brain"), " AlphaFold Structure"), tags$br(),
                      tags$small("Full protein coverage (computational prediction). Enable 'Color by pLDDT' to see confidence.")))
    }
    
    identity_str <- if (!is.null(rv$identity) && !is.na(rv$identity)) paste0(round(rv$identity * 100, 1), "%") else "N/A"
    
    tags$div(class = "alert alert-success", style = "margin-bottom: 10px;",
             tags$strong(icon("cube"), paste(" PDB:", toupper(rv$pdb_id), "(chain", rv$chain_id, ")")), tags$br(),
             tags$small(paste("UniProt:", rv$unp_start, "-", rv$unp_end, "| Identity:", identity_str)))
  })
  
  # Mutation warning panel
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
      from_aa <- if (!is.null(m$mutation_details$from)) m$mutation_details$from else "?"
      to_aa <- if (!is.null(m$mutation_details$to)) m$mutation_details$to else "?"
      resnum <- if (!is.null(m$residue_number)) m$residue_number else "?"
      paste0(from_aa, resnum, to_aa)
    })
    
    tags$div(class = "mutation-warning",
             tags$strong(icon("exclamation-triangle"), " ENGINEERED MUTATIONS DETECTED (Chain ", rv$chain_id, ")"),
             tags$br(), tags$br(),
             tags$span(style = "font-size: 14px; color: #856404; font-weight: bold;", paste(mut_texts, collapse = ", ")),
             tags$br(), tags$br(),
             tags$span(style = "font-size: 11px;",
                       "These residues were ", tags$strong("mutated by crystallographers"), " to improve stability or crystallization. ",
                       "This structure may ", tags$strong("NOT represent wild-type"), ". ",
                       "Mutations shown in ", tags$span(style = "color: orange; font-weight: bold;", "ORANGE"), " in the 3D viewer.")
    )
  })
  
  # Assembly info panel
  output$assembly_info <- renderUI({
    if (is.null(rv$pdb_summary) || rv$is_alphafold) return(NULL)
    
    summ <- rv$pdb_summary
    assemblies <- summ$assemblies
    n_assemblies <- if (!is.null(assemblies)) length(assemblies) else 0
    title <- summ$title %||% "No title"
    
    tags$div(class = "assembly-info",
             tags$strong(icon("puzzle-piece"), " Structure Information"),
             tags$br(), tags$br(),
             tags$strong("Title: "), tags$em(title),
             tags$br(),
             if (n_assemblies > 0) {
               tagList(
                 tags$strong("Biological Assemblies: "), n_assemblies,
                 tags$br(),
                 tags$small("The biological assembly represents the functional oligomeric state. ",
                            "Variants at interfaces may disrupt protein-protein interactions.")
               )
             }
    )
  })
  
  # NGL Viewer
  output$ngl_view <- NGLVieweR::renderNGLVieweR({
    req(rv$pdb_text)
    
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
            sele <- paste0(dom$start, "-", dom$end, ":A")
            color <- DOMAIN_COLORS[(i %% length(DOMAIN_COLORS)) + 1]
            ngl <- ngl %>% addRepresentation("cartoon", param = list(sele = sele, color = color, opacity = 1.0))
          }
        }
      } else {
        ngl <- ngl %>% addRepresentation("cartoon", param = list(sele = "protein", color = "spectrum"))
      }
    } else {
      # PDB structure
      if (input$show_domains && !is.null(rv$domains) && nrow(rv$domains) > 0) {
        ngl <- ngl %>%
          addRepresentation("cartoon", param = list(sele = paste0(":", rv$chain_id, " and protein"), color = "lightgray", opacity = 0.5)) %>%
          addRepresentation("cartoon", param = list(sele = paste0("not :", rv$chain_id, " and protein"), color = "lightgray", opacity = 0.15))
        
        for (i in seq_len(min(nrow(rv$domains), 10))) {
          dom <- rv$domains[i, ]
          if (!is.na(dom$start) && !is.na(dom$end)) {
            if (dom$end >= rv$unp_start && dom$start <= rv$unp_end) {
              dom_pdb_start <- rv$pdb_start + max(0, dom$start - rv$unp_start)
              dom_pdb_end <- rv$pdb_start + min(dom$end - rv$unp_start, rv$unp_end - rv$unp_start)
              sele <- paste0(dom_pdb_start, "-", dom_pdb_end, ":", rv$chain_id)
              color <- DOMAIN_COLORS[(i %% length(DOMAIN_COLORS)) + 1]
              ngl <- ngl %>% addRepresentation("cartoon", param = list(sele = sele, color = color, opacity = 1.0))
            }
          }
        }
      } else {
        ngl <- ngl %>%
          addRepresentation("cartoon", param = list(sele = paste0(":", rv$chain_id, " and protein"), color = "steelblue")) %>%
          addRepresentation("cartoon", param = list(sele = paste0("not :", rv$chain_id, " and protein"), color = "lightgray", opacity = 0.2))
      }
      
      # Engineered mutations in ORANGE
      if (!is.null(rv$mutations) && length(rv$mutations) > 0) {
        for (mut in rv$mutations) {
          if (!is.null(mut$chain_id) && mut$chain_id == rv$chain_id && !is.null(mut$residue_number)) {
            sele <- paste0(mut$residue_number, ":", rv$chain_id)
            ngl <- ngl %>% addRepresentation("ball+stick", param = list(sele = sele, color = "orange", radius = 0.3))
          }
        }
      }
    }
    
    # Show ligands if checkbox is checked
    if (input$show_ligands) {
      ngl <- ngl %>% addRepresentation("ball+stick", 
                                       param = list(name = "ligands_init", sele = "ligand", colorScheme = "element"))
    }
    
    # Start with rocking ON by default
    ngl <- ngl %>% setRock(TRUE)
    
    ngl
  })
  
  # Spinning control
  observeEvent(input$spinning, {
    if (input$spinning) {
      NGLVieweR_proxy("ngl_view") %>% updateSpin(TRUE)
      if (input$rocking) updateCheckboxInput(session, "rocking", value = FALSE)
    } else {
      NGLVieweR_proxy("ngl_view") %>% updateSpin(FALSE)
    }
  }, ignoreInit = TRUE)
  
  # Rocking control - when unchecked, STOP rocking
  observeEvent(input$rocking, {
    if (input$rocking) {
      NGLVieweR_proxy("ngl_view") %>% updateRock(TRUE)
      if (input$spinning) updateCheckboxInput(session, "spinning", value = FALSE)
    } else {
      NGLVieweR_proxy("ngl_view") %>% updateRock(FALSE)
    }
  }, ignoreInit = TRUE)
  
  # Ramachandran plot
  output$rama_plot <- renderPlotly({
    tor <- rv$torsion_data
    if (is.null(tor)) return(NULL)
    
    highlight <- rv$selected_variant_pos
    if (!is.null(rv$selected_rama_residue)) highlight <- rv$selected_rama_residue
    
    title <- paste("Ramachandran -", toupper(rv$pdb_id %||% ""))
    create_ramachandran_plotly(tor, highlight_resno = highlight, title = title)
  })
  
  observeEvent(event_data("plotly_click", source = "rama_plot"), {
    click <- event_data("plotly_click", source = "rama_plot")
    if (!is.null(click) && !is.null(click$key)) {
      rv$selected_rama_residue <- as.integer(click$key)
    }
  })
  
  output$rama_click_info <- renderPrint({
    if (is.null(rv$selected_rama_residue)) { cat("Click a point to see details"); return() }
    
    tor <- rv$torsion_data
    if (is.null(tor)) return()
    
    res <- tor[tor$resno == rv$selected_rama_residue, ]
    if (nrow(res) == 0) { cat("Click a point to see details"); return() }
    
    cat("Selected:", res$resid[1], res$resno[1], "\n")
    cat("Chain:", res$chain[1], "\n")
    cat("φ:", round(res$phi[1], 1), "° | ψ:", round(res$psi[1], 1), "°\n")
    cat("Category:", res$category[1])
  })
  
  outputOptions(output, "ngl_view", suspendWhenHidden = FALSE)
}

shinyApp(ui = ui, server = server)