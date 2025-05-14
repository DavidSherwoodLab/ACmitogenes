# ============================================================
# Anchor Cell Mitochondrial Gene Transcriptome Explorer
# Author: Jake Leyhr
# Contact: jake.leyhr@duke.edu
# Created: 2025-05-13
# Description:
#     A Shiny app for interactive exploration of mitochondrial
#     gene expression and orthology in the C. elegans 
#     anchor cell, including GSEA pathway visualization.
#
# License: MIT / CC-BY
# ============================================================

# Load required libraries
library(shiny)
library(shinyjs)
library(data.table)
library(tidyverse)
library(DT)
library(scales)
library(plotly)
library(forcats)
library(heatmaply)

# Load main expression + orthology data
merged_data <- fread("merged_orthologs_mitopathways_transcriptome.csv")

# Load and structure GSEA data — split into nested pathway levels
gsea_data <- fread("gsea_all_results.csv") %>%
  mutate(
    Significance = ifelse(FDR.q.val < 0.05, "Significant", "Not Significant"),
    Level1 = sapply(strsplit(NAME, "_>_"), `[`, 1),
    Level2 = sapply(strsplit(NAME, "_>_"), function(x) if(length(x) >= 2) x[2] else NA),
    Level3 = sapply(strsplit(NAME, "_>_"), function(x) if(length(x) >= 3) x[3] else NA)
  )

# Shiny UI
ui <- fluidPage(
  useShinyjs(),
  tagList(
    titlePanel("Anchor Cell Mitochondrial Gene Transcriptome Explorer"),
    tags$div(
      HTML("Interactively explore expression changes of mitochondrial pathway orthologs in the <i>C. elegans</i> anchor cell."),
      tags$br(),
      "Read the associated publication: ",
      tags$a(
        href = "https://doi.org/10.1101/2025.05.02.651978",
        "Kenny-Ganzert et al. (2025) Specialized high-capacity mitochondria fuel cell invasion. bioRxiv",
        target = "_blank"
      )
    )
  ),
  
  
  tabsetPanel(
    tabPanel("Gene & Pathway Search",
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 helpText("Search by gene name or pathway."),
                 textInput("human_gene", "Enter Human Gene Names (comma or space separated):", value = ""),
                 textInput("worm_gene", "Enter Worm Gene Names (comma or space separated):", value = ""),
                 selectizeInput("search_pathway", "Select Pathway (choose from dropdown or start typing):", choices = NULL, multiple = TRUE),
                 actionButton("search", "Search"),
                 hr(),
                 downloadButton("download_table", "Download Genes Table")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Volcano Plot", plotlyOutput("volcano_plot"), height = "600px", width = "800px"),
                   tabPanel("Bar Plot", plotlyOutput("bar_chart", height = "600px"),
                            checkboxInput("log_scale", "Log scale y-axis", value = FALSE),
                            checkboxInput("order_by_ac", "Order by AC expression (high to low)", value = FALSE)
                   ),
                   tabPanel("Genes Table", DTOutput("filtered_table"))
                 )
               )
             )
    ),
    tabPanel("Explore GSEA Pathways",
             sidebarLayout(
               sidebarPanel(
                 width = 12,
                 fluidRow(
                   column(6,
                          h4("Navigate GSEA Pathways"),
                          helpText("Click a bar to enter into sub-pathways.")
                   ),
                   column(6, align = "right",
                          div(style = "margin-top: 15px;",
                              actionButton("gsea_back", "⬅ Go Back"),
                              actionButton("gsea_reset", "⏮ Reset")
                          ))),
                 br(),
                 plotlyOutput("gsea_nested_barplot", height = "500px", width = "100%")
               ),
               mainPanel(
                 width = 12,
                 actionButton("explore_genes", "Explore Genes in Selected Pathway"),
                 downloadButton("download_table2", "Download Genes Table"),
                 tabsetPanel(
                   tabPanel("Volcano Plot", plotlyOutput("selected_volcano", height = "600px")),
                   tabPanel("Bar Plot", plotlyOutput("selected_barplot", height = "600px"),
                            checkboxInput("log_scale", "Log scale y-axis", value = FALSE),
                            checkboxInput("order_by_ac", "Order by AC expression (high to low)", value = FALSE)),
                   tabPanel("Heatmap", plotlyOutput("selected_heatmap", height = "600px"),
                            checkboxInput("filter_significant", "Only show significantly differentially expressed genes (padj < 0.05)", value = FALSE),
                            checkboxInput("cluster_by_path", "Cluster genes by pathway", value = FALSE),
                   ),
                   tabPanel("Genes Table", DTOutput("selected_pathway_table"))
                 )
               )
             )
    ),
    tabPanel("Gene Homologies",
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 helpText("Search for orthologs by gene name."),
                 textInput("human_gene", "Enter Human Gene Names (comma or space separated):", value = ""),
                 textInput("worm_gene", "Enter Worm Gene Names (comma or space separated):", value = ""),
                 actionButton("search", "Search"),
                 textOutput("error_message")
               ),
               mainPanel(
                 h4("Gene Homology Matches Across Databases", style = "margin-bottom: 70px;"),
                 
                 # Wrap in unique ID for scoped CSS
                 div(id = "homology_table_wrapper",
                     DTOutput("homologies_table")
                 ),
                 
                 # Scoped CSS for this table only
                 tags$head(
                   tags$style(HTML("
                    /* Scoped to only this table */
                    #homology_table_wrapper .dataTable thead th {
                      transform: rotate(-45deg);
                      transform-origin: bottom left;
                      white-space: nowrap !important;
                      text-align: left !important;
                      vertical-align: bottom !important;
                      padding: 2px 4px !important;
                      font-size: 12px !important;
                      min-width: 30px !important;
                      max-width: 60px !important;
                    }
                  
                    #homology_table_wrapper .dataTable td:first-child,
                    #homology_table_wrapper .dataTable th:first-child,
                    #homology_table_wrapper .dataTable td:nth-child(2),
                    #homology_table_wrapper .dataTable th:nth-child(2) {
                      min-width: 200px !important;
                      max-width: 400px !important;
                      padding-left: 10px !important;
                      padding-right: 10px !important;
                      white-space: normal !important;
                    }
                  
                    #homology_table_wrapper .dataTable td {
                      padding: 2px 4px !important;
                      font-size: 14px !important;
                    }
                  
                    #homology_table_wrapper .dataTables_wrapper {
                      font-size: 12px;
                    }
                  
                    #homology_table_wrapper table.dataTable {
                      table-layout: fixed !important;
                    }
                  
                    /* Optional: Alternating color for algorithm columns */
                    #homology_table_wrapper table.dataTable thead th:nth-child(n+3):nth-child(odd),
                    #homology_table_wrapper table.dataTable tbody td:nth-child(n+3):nth-child(odd) {
                      background-color: #f0f0f0;
                    }
                  
                    #homology_table_wrapper table.dataTable thead th:nth-child(n+3):nth-child(even),
                    #homology_table_wrapper table.dataTable tbody td:nth-child(n+3):nth-child(even) {
                      background-color: #ffffff;
                    }
                  "))
                 )
               )
             )
    )
  )
)


# Server logic
server <- function(input, output, session) {
  updateSelectizeInput(session, "search_pathway", choices = unique(merged_data$full_pathway), server = TRUE)
  
  # =======================
  # Gene & Pathway Search Tab
  # =======================
  # Event-triggered reactive filtering — based on human/worm gene or pathway input
  filtered_data <- eventReactive(input$search, {
    # Parse input into vectors of search terms
    worm_input <- trimws(unlist(strsplit(input$worm_gene, "[[:space:],]+")))
    human_input <- trimws(unlist(strsplit(input$human_gene, "[[:space:],]+")))
    path_input <- input$search_pathway
    
    # Ensure only one input field is active
    if ((length(worm_input) > 0) + (length(human_input) > 0) + (length(path_input) > 0) > 1) {
      validate(need(FALSE, "Please provide input in only one field."))
    }
    if (length(worm_input) == 0 && length(human_input) == 0 && length(path_input) == 0) {
      validate(need(FALSE, "Please enter gene names or select a pathway."))
    }
    
    # Apply filtering logic
    filtered <- if (length(worm_input) > 0) {
      merged_data %>%
        filter(
          tolower(worm_gene) %in% tolower(worm_input) |
            tolower(`Updated Gene Name`) %in% tolower(worm_input) |
            tolower(worm_ID) %in% tolower(worm_input)
        )
    } else if (length(human_input) > 0) {
      merged_data %>% filter(tolower(human_gene) %in% tolower(human_input))
    } else {
      merged_data %>% filter(full_pathway %in% path_input | Pathway_Level %in% path_input)
    }
    
    # Drop duplicated genes and remove unused columns
    filtered %>%
      distinct(`Updated Gene Name`, .keep_all = TRUE) %>%
      select(-Pathway_Level, -geneName, -`Updated Gene Name`)
  })
  
  # Make volcano plot
  output$volcano_plot <- renderPlotly({
    df <- filtered_data()
    
    # Calculate mean expression across replicates
    df <- df %>% mutate(
      mean_AC = rowMeans(select(., AC1, AC2, AC3)),
      mean_WB = rowMeans(select(., WB1, WB2, WB3))
    )
    
    # Classify genes by log2FC and padj
    df <- df %>% mutate(
      diffexpressed = case_when(
        `log2 fold change` > 1 & padj < 0.05 ~ "UP",
        `log2 fold change` < -1 & padj < 0.05 ~ "DOWN",
        TRUE ~ "NO"
      )
    )
    
    # Build plot
    p <- ggplot(df, aes(
      x = `log2 fold change`,
      y = -log10(padj),
      color = diffexpressed,
      text = paste0(
        "Gene: ", worm_gene,
        "<br>log2FC: ", round(`log2 fold change`, 2),
        "<br>padj: ", signif(padj, 3),
        "<br>Mean AC Reads: ", round(mean_AC),
        "<br>Mean WB Reads: ", round(mean_WB)
      )
    )) +
      geom_point(aes(size = mean_AC)) +
      scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NO" = "grey")) +
      guides(size = "none") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_classic() +
      labs(title = "AC Volcano Plot", x = "log2 Fold Change", y = "-log10(padj)")
    
    ggplotly(p, tooltip = "text")
  })
  
  # Make bar plot
  output$bar_chart <- renderPlotly({
    df <- filtered_data()
    df <- df %>% mutate(
      mean_AC = rowMeans(select(., AC1, AC2, AC3)),
      mean_WB = rowMeans(select(., WB1, WB2, WB3))
    )
    
    if (input$order_by_ac) {
      df <- df %>% mutate(worm_gene = fct_reorder(worm_gene, mean_AC, .desc = TRUE))
    } else {
      df <- df %>% mutate(worm_gene = factor(worm_gene, levels = sort(unique(worm_gene))))
    }
    
    df_long <- df %>%
      pivot_longer(cols = c(mean_AC, mean_WB), names_to = "Tissue", values_to = "Expression") %>%
      mutate(Tissue = recode(Tissue, mean_AC = "AC", mean_WB = "WB"))
    
    p <- ggplot(df_long, aes(
      x = worm_gene,
      y = Expression,
      fill = Tissue,
      text = paste0("Gene: ", worm_gene, "<br>Tissue: ", Tissue, "<br>Mean Reads: ", round(Expression, 2))
    )) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("AC" = "red", "WB" = "lightgreen")) +
      theme_minimal() +
      labs(title = "AC Expression", x = "Gene", y = "Expression (reads)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (input$log_scale) {
      p <- p + scale_y_log10(limits = c(1, NA), expand = c(0, 0))
    } else {
      p <- p + scale_y_continuous(expand = c(0, 0))
    }
    
    ggplotly(p, tooltip = "text")
  })
  
  # Make Genes Table  
  output$filtered_table <- renderDT({
    datatable(filtered_data())
  })
  
  # Function to download Genes Table
  output$download_table <- downloadHandler(
    filename = function() "filtered_genes.csv",
    content = function(file) write.csv(filtered_data(), file, row.names = FALSE)
  )
  
  
  
  # =======================
  # Explore GSEA Pathways Tab
  # =======================
  # Stores current level of pathway nesting (Level1 → Level2 → Level3) and selected pathways at each level
  gsea_path <- reactiveValues(level = "Level1", selections = list())
  
  # Holds currently selected terminal pathway (i.e., one without sub-pathways)
  selected_terminal <- reactiveVal(NULL)
  
  # Handle clicks on the GSEA bar plot
  observeEvent(event_data("plotly_click", source = "gsea_plot"), {
    d <- event_data("plotly_click", source = "gsea_plot")
    if (!is.null(d)) {
      bar_id <- d$customdata
      if (is.null(bar_id) || (length(bar_id) == 1 && is.na(bar_id))) {
        return()
      }
      
      # Remove trailing index suffix from clicked bar label (e.g., "Pathway_3" → "Pathway")
      sel <- sub("_\\d+$", "", bar_id)
      if (is.na(sel) || sel == "NA") {
        sel <- "__NA__"
      }
      
      lvl <- gsea_path$level
      
      # Check if this selection has sub-pathways at the next level
      has_children <- if (lvl == "Level1") {
        any(gsea_data$Level1 == sel & !is.na(gsea_data$Level2) & gsea_data$Level2 != "")
      } else if (lvl == "Level2") {
        any(gsea_data$Level2 == sel & !is.na(gsea_data$Level3) & gsea_data$Level3 != "")
      } else {
        FALSE
      }
      
      # If the selection has children, move down one level
      if (has_children) {
        if (lvl == "Level1") {
          gsea_path$selections$Level1 <- sel
          gsea_path$level <- "Level2"
        } else if (lvl == "Level2") {
          gsea_path$selections$Level2 <- sel
          gsea_path$level <- "Level3"
        }
        selected_terminal(NULL) # Reset terminal selection
      } else {
        # Toggle selection of the terminal pathway (click again to deselect)
        if (identical(selected_terminal(), sel)) {
          selected_terminal(NULL)
        } else {
          selected_terminal(sel)
        }
      }
    }
  })
  
  # Handle "⬅ Go Back" — move up one level in the hierarchy
  observeEvent(input$gsea_back, {
    if (gsea_path$level == "Level3") {
      gsea_path$level <- "Level2"
      gsea_path$selections$Level2 <- NULL
    } else if (gsea_path$level == "Level2") {
      gsea_path$level <- "Level1"
      gsea_path$selections$Level1 <- NULL
    }
    # Reset the terminal
    selected_terminal(NULL)
  })
  
  # Handle "⏮ Reset" — go all the way back to Level1
  observeEvent(input$gsea_reset, {
    gsea_path$level <- "Level1"
    gsea_path$selections <- list()
    selected_terminal(NULL)
  })
  
  # Render the nested GSEA barplot for the current level
  output$gsea_nested_barplot <- renderPlotly({
    level <- gsea_path$level
    sel1 <- gsea_path$selections$Level1
    sel2 <- gsea_path$selections$Level2
    
    # Get relevant data for the current nesting level
    if (level == "Level1") {
      df <- gsea_data %>%
        filter(!is.na(Level1) & (is.na(Level2) | Level2 == "")) %>%
        mutate(
          n_next = map_int(Level1, ~ n_distinct(gsea_data$Level2[
            !is.na(gsea_data$Level1) & gsea_data$Level1 == .x &
              !is.na(gsea_data$Level2) & gsea_data$Level2 != ""
          ]))
        )
      col <- "Level1"
        } else if (level == "Level2") {
          df <- gsea_data %>%
            filter(Level1 == sel1 & (is.na(Level3) | Level3 == "")) %>%
            mutate(
              n_next = map_int(Level2, function(parent) {
                if (is.na(parent)) return(0)
                n_distinct(gsea_data$Level3[
                  !is.na(gsea_data$Level2) & gsea_data$Level2 == parent &
                    !is.na(gsea_data$Level3) & gsea_data$Level3 != ""
                ])
              })
            )
          col <- "Level2"
        } else {
          df <- gsea_data %>%
            filter(Level1 == sel1, Level2 == sel2, !is.na(Level3) & Level3 != "") %>%
            mutate(n_next = 0)
          col <- "Level3"
    }
    
    # Format labels and assign colors
    df$label <- paste0(df[[col]], " (", df$n_next, ")")
    df$fill_color <- ifelse(df$NES < 0, "steelblue", "lightblue")
    df$label <- paste0(df[[col]], " (", df$n_next, ")")
    
    if (level == "Level3" && !is.null(sel1) && !is.null(sel2)) {
      path_title <- paste(sel1, ">", sel2)
    } else if (level == "Level2" && !is.null(sel1)) {
      path_title <- sel1
    } else {
      path_title <- ""
    }
    
    # Ensure the column is character
    df[[col]] <- as.character(df[[col]])
    # Trim whitespace
    df[[col]] <- trimws(df[[col]])
    df[[col]] <- ifelse(is.na(df[[col]]), "__NA__", df[[col]])
    
    # Create hover text for each bar 
    df$hover_text <- paste0(
      df[[col]], "<br>",
      "AC NES: ", signif(df$NES, 3), "<br>",
      "Sub-pathways: ", df$n_next, "<br>",
      "FDR q-value: ", signif(df$FDR.q.val, 3)
    )
    
    # Assign a unique ID to each bar
    df$bar_id <- paste0(df[[col]], "_", seq_len(nrow(df)))
    
    # Sort pathways by enrichment score
    df <- df %>% arrange(NES)
    
    # Build display label for each bar
    df$label <- paste0(df[[col]], " (", df$n_next, ")")
    
    # Convert label to factor with levels in the order they appear
    df$label <- factor(df$label, levels = df$label)
    
    # Highlight the selected terminal bar
    selected <- selected_terminal()
    df$line_width <- ifelse(sub("_\\d+$", "", df$bar_id) == selected, 5, 0)
    df$line_color <- ifelse(df$line_width > 0, "red", "transparent")
    
    # Plot the pathways bar chart
    plot_ly(
      data = df,
      x = ~NES,
      y = ~label,
      type = 'bar',
      orientation = 'h',
      marker = list(
        color = df$fill_color,
        line = list(color = df$line_color, width = df$line_width)
      ),
      customdata = ~bar_id,
      source = "gsea_plot",
      text = I(df$hover_text),  # force plotly to use this literal text
      textposition = "none",
      hovertemplate = "%{text}<extra></extra>"
    ) %>%
      plotly::event_register("plotly_click") %>%
      layout(
        title = list(
          text = paste("Pathway:", path_title),
          y = 0.95  # move title down slightly to leave more space for the menu
        ),
        margin = list(t = 80),  # add top padding
        yaxis = list(title = "", fixedrange = TRUE),
        xaxis = list(title = "AC Normalized Enrichment Score", fixedrange = TRUE, 
            range = c(
              min(min(df$NES, na.rm = TRUE), -1),
              max(max(df$NES, na.rm = TRUE), 1)
            )),
        dragmode = FALSE  # This disables drag-to-zoom and ensures click is captured immediately
        )  %>%
      config(
        displayModeBar = TRUE,
        modeBarButtonsToRemove = c(
          "zoom2d", "pan2d", "select2d", "lasso2d",
          "zoomIn2d", "zoomOut2d", "autoScale2d", "resetScale2d"
        ),
        scrollZoom = FALSE
      )
  })
  
  # Reactive expression to return the most specific selected pathway
  # Priority order: Level3 > Level2 > Level1
  selected_path <- reactive({
    sel1 <- gsea_path$selections$Level1
    sel2 <- gsea_path$selections$Level2
    sel3 <- selected_terminal()
    if (!is.null(sel3)) {
      sel3
    } else if (!is.null(sel2)) {
      sel2
    } else if (!is.null(sel1)) {
      sel1
    } else {
      NULL
    }
  })
  
  # Stores the most recently confirmed pathway selection (via button click)
  confirmed_path <- reactiveVal(NULL)
  
  # When "Explore Genes" button is clicked, lock in the current pathway selection
  observeEvent(input$explore_genes, {
    confirmed_path(selected_path())
  })
  
  # Returns a subset of genes from merged_data based on the selected pathway
  selected_genes <- eventReactive(input$explore_genes, {
    path <- selected_path()
    # cat("DEBUG: selected_path =", path, "\n")
    # if (is.null(path)) {
    #   cat("No valid selection — returning empty tibble\n")
    #   return(tibble())
    # }
    
    if (path == "__NA__") {
      sel1 <- gsea_path$selections$Level1
      # cat("DEBUG: Handling '__NA__' bar click\n")
      # cat("DEBUG: Level1 context =", sel1, "\n")
      
      # Normalize for comparison (lowercase, trim, replace underscores)
      norm_sel1 <- trimws(tolower(gsub("_", " ", sel1)))
      
      # Filter genes matching this Level1 label and ensure it's not a nested pathway
      result <- merged_data %>%
        filter(
          !is.na(Pathway_Level) &
            trimws(tolower(Pathway_Level)) == norm_sel1 & 
            # exclude entries with deeper nesting
            !grepl(">", full_pathway, fixed = TRUE)
        )
      
      # cat("DEBUG: Rows matched =", nrow(result), "\n")
      # cat("DEBUG: Sample matched Pathway_Level:\n")
      # print(unique(result$Pathway_Level))
      return(result)
    }
    
    # General case: search for the selected pathway in the Pathway_Level column
    search_string <- tolower(gsub("_", " ", path))
    # cat("DEBUG: Searching case-insensitive with word boundary for:", search_string, "\n")
    result <- merged_data %>%
      filter(grepl(paste0("\\b", search_string, "\\b"), tolower(Pathway_Level), ignore.case = TRUE))
    # cat("DEBUG: Rows matched =", nrow(result), "\n")
    # cat("DEBUG: Sample matched Pathway_Level:\n")
    # print(head(result$Pathway_Level))
    return(result)
  })

  # Make the volcano plot
  output$selected_volcano <- renderPlotly({
    df <- selected_genes()
    df <- df %>% mutate(
      mean_AC = rowMeans(select(., AC1, AC2, AC3)),
      mean_WB = rowMeans(select(., WB1, WB2, WB3))
    )
    req(nrow(df) > 0)
    df <- df %>% mutate(
      diffexpressed = case_when(
        `log2 fold change` > 1 & padj < 0.05 ~ "UP",
        `log2 fold change` < -1 & padj < 0.05 ~ "DOWN",
        TRUE ~ "NO"
      )
    )
    
    p <- ggplot(df, aes(
      x = `log2 fold change`,
      y = -log10(padj),
      color = diffexpressed,
      text = paste0(
        "Gene: ", worm_gene,
        "<br>log2FC: ", round(`log2 fold change`, 2),
        "<br>padj: ", signif(padj, 3),
        "<br>Mean AC Reads: ", round(mean_AC),
        "<br>Mean WB Reads: ", round(mean_WB)
      )
    )) +
      geom_point(aes(size = mean_AC)) +
      scale_color_manual(
        name = "Differential Expression",
        values = c("UP" = "red", "DOWN" = "blue", "NO" = "grey")
      ) +
      guides(size = "none") +  # hides the size legend
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_classic() +
      labs(title = paste(confirmed_path(), "AC Volcano Plot"),
           x = "log2 Fold Change", y = "-log10(padj)") 
  ggplotly(p, tooltip = "text")
  })

  
  # Make the barplot
  output$selected_barplot <- renderPlotly({
    df <- selected_genes()
    df <- df %>% mutate(
      mean_AC = rowMeans(select(., AC1, AC2, AC3)),
      mean_WB = rowMeans(select(., WB1, WB2, WB3))
    )
    
    # Order genes based on AC expression if checkbox is selected
    if (input$order_by_ac) {
      df <- df %>% mutate(worm_gene = fct_reorder(worm_gene, mean_AC, .desc = TRUE))
    } else {
      df <- df %>% mutate(worm_gene = factor(worm_gene, levels = sort(unique(worm_gene))))
    }
    
    # Reshape data to long format for bar plot (one row per gene × tissue)
    df_long <- df %>%
      pivot_longer(cols = c(mean_AC, mean_WB), names_to = "Tissue", values_to = "Expression") %>%
      mutate(Tissue = recode(Tissue, mean_AC = "AC", mean_WB = "WB"))
    
    p <- ggplot(df_long, aes(
      x = worm_gene,
      y = Expression,
      fill = Tissue,
      text = paste0("Gene: ", worm_gene, "<br>Tissue: ", Tissue, "<br>Mean Reads: ", round(Expression, 2))
    )) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("AC" = "red", "WB" = "lightgreen")) +
      theme_minimal() +
      labs(title = "AC Expression", x = "Gene", y = "Expression (reads)") +
      labs(title = paste(confirmed_path(), "AC Expression"), y = "Expression (reads)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Optionally apply log scale to y-axis
    if (input$log_scale) {
      p <- p + scale_y_log10(limits = c(1, NA), expand = c(0, 0))
    } else {
      p <- p + scale_y_continuous(expand = c(0, 0))
    }
    
    ggplotly(p, tooltip = "text")
  })
  
  # Make Genes Table
  output$selected_pathway_table <- renderDT({
    datatable(selected_genes())
  })
  
  # Function to download Genes Table
  output$download_table2 <- downloadHandler(
    filename = function() "filtered_genes.csv",
    content = function(file) write.csv(selected_genes(), file, row.names = FALSE)
  )
  
  # Make heatmap
  output$selected_heatmap <- renderPlotly({
    df <- selected_genes()
  
    # Optionally filter to only significantly differentially expressed genes
    if (input$filter_significant) {
      df <- df %>% filter(!is.na(padj) & padj < 0.05)
    }
    # Exit if no rows remain after filtering
    req(nrow(df) > 0)
    
    # Average replicates for each gene × tissue combination
    df <- df %>%
      group_by(worm_gene) %>%
      summarise(
        AC1 = mean(AC1, na.rm = TRUE),
        AC2 = mean(AC2, na.rm = TRUE),
        AC3 = mean(AC3, na.rm = TRUE),
        WB1 = mean(WB1, na.rm = TRUE),
        WB2 = mean(WB2, na.rm = TRUE),
        WB3 = mean(WB3, na.rm = TRUE),
        Pathway_Level = last(Pathway_Level),
        .groups = "drop"
      ) %>% 
    arrange(Pathway_Level) 
    
    # Apply log10 transform (add 1 to avoid log(0))
    df <- df %>%
      mutate(across(c(AC1, AC2, AC3, WB1, WB2, WB3), ~log10(.x + 1)))
    
    # Reshape into matrix format for heatmaply
    expr_matrix <- df %>%
      pivot_longer(cols = c(AC1, AC2, AC3, WB1, WB2, WB3), names_to = "Sample", values_to = "Expression") %>%
      pivot_wider(names_from = Sample, values_from = Expression, values_fn = list(Expression = mean)) %>%
      select(worm_gene, AC1, AC2, AC3, WB1, WB2, WB3) %>%  # force correct column order
      mutate(across(-worm_gene, ~ as.numeric(trimws(.)))) %>%  # convert all to numeric        
      as.data.frame() %>%
      column_to_rownames("worm_gene") %>%
      as.matrix()
    
    # Set up pathway-level row annotations
    row_ann <- df %>%
      select(worm_gene, Pathway_Level) %>%
      distinct() %>%
      as.data.frame() %>%
      column_to_rownames("worm_gene")
    
    # Generate the heatmap with clustering and annotations   
    heatmaply::heatmaply(
      expr_matrix,
      Rowv = if (input$cluster_by_path) NA else TRUE,
      Colv = NA,
      row_side_colors = row_ann,
      scale = "row",
      colors = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
      colorbar = FALSE,
      show_colorbar = FALSE,
      main = paste("Expression Heatmap:", confirmed_path()),
      xlab = "Sample", ylab = "Worm Gene",
    )
  })
  
  
  # =======================
  # Gene Homologies Tab
  # =======================
  # Filter orthology table based on user input (triggered by "Search" button)
  filtered_homologies <- eventReactive(input$search, {
    # Parse input: split on spaces or commas, trim whitespace
    worm_input <- trimws(unlist(strsplit(input$worm_gene, "[[:space:],]+")))
    human_input <- trimws(unlist(strsplit(input$human_gene, "[[:space:],]+")))
    
    # Ensure only one input field is active
    if ((length(worm_input) > 0) + (length(human_input) > 0) > 1) {
      validate(need(FALSE, "Please provide input in only one field."))
    }
    if (length(worm_input) == 0 && length(human_input) == 0) {
      validate(need(FALSE, "Please enter a gene name to search."))
    }
    
    # Filter merged_data based on user input
    df <- if (length(worm_input) > 0) {
      merged_data %>%
        filter(
          tolower(worm_gene) %in% tolower(worm_input) |
            tolower(worm_ID) %in% tolower(worm_input)
        )
    } else {
      merged_data %>%
        filter(tolower(human_gene) %in% tolower(human_input))
    }
    
    # List of orthology databases to display as individual columns
    known_databases <- c(
      "Ensembl Compara", "Hieranoid", "Homologene", "InParanoid", "OMA",
      "OrthoFinder", "OrthoInspector", "OrthoMCL", "PANTHER",
      "PhylomeDB", "SonicParanoid"
    )
    
    
    for (alg in known_databases) {
      df[[alg]] <- ifelse(str_detect(df$all_databases, fixed(alg)), "✓", "")
    }
    
    colnames(df)[colnames(df) == "human_gene"] <- "Human"
    colnames(df)[colnames(df) == "worm_gene"] <- "Worm"
    
    
    df %>%
      select(Human, Worm, all_of(known_databases)) %>%
      distinct()
  })
  
  
  
  output$homologies_table <- DT::renderDataTable({
    DT::datatable(
      filtered_homologies(),
      options = list(
        pageLength = 25,
        dom = 't',  # table only
        ordering = FALSE,
        columnDefs = list(
          list(width = '100px', targets = 0),  # Human column
          list(width = '100px', targets = 1)   # Worm column
        )
      ),
      rownames = FALSE
    )
  })
  

}

# Run app
shinyApp(ui = ui, server = server)

