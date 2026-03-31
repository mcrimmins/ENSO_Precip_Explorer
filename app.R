library(shiny)
library(shinycssloaders)
library(tidyverse)
library(sf)
library(leaflet)
library(rvest)
library(janitor)
library(ggrepel)
library(glue)
library(scales)
library(plotly)
library(data.table)

options(shiny.maxRequestSize = 100 * 1024^2)

# -----------------------------
# Configuration
# -----------------------------
DATA_DIR <- "data"
CACHE_DIR <- "cache"
SHAPEFILE_NAME <- "GIS.OFFICIAL_CLIM_DIVISIONS.shp"
CLIMDIV_BASE_URL <- "https://www.ncei.noaa.gov/pub/data/cirs/climdiv"
RONI_URL <- "https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt"
PCP_PREFIX <- "climdiv-pcpndv-v1.0.0"
RONI_CACHE_FILE <- "RONI.ascii.txt"

SEASON_MONTHS <- list(
  DJF = c(12, 1, 2),
  JFM = c(1, 2, 3),
  FMA = c(2, 3, 4),
  MAM = c(3, 4, 5),
  AMJ = c(4, 5, 6),
  MJJ = c(5, 6, 7),
  JJA = c(6, 7, 8),
  JAS = c(7, 8, 9),
  ASO = c(8, 9, 10),
  SON = c(9, 10, 11),
  OND = c(10, 11, 12),
  NDJ = c(11, 12, 1)
)

ENSO_LEVELS <- c("La Nina", "Neutral", "El Nino")
ENSO_COLORS <- c(
  "La Nina" = "#1f4ed8",
  "Neutral" = "#20c933",
  "El Nino" = "#ff1f1f"
)

# -----------------------------
# Utilities
# -----------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

elapsed_secs <- function(t0) {
  round(as.numeric(Sys.time() - t0, units = "secs"), 2)
}

ensure_cache_dir <- function(cache_dir = CACHE_DIR) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  invisible(cache_dir)
}

get_latest_procdate <- function(base_url = CLIMDIV_BASE_URL) {
  procdate_url <- paste0(base_url, "/procdate.txt")
  procdate <- tryCatch(
    readr::read_lines(procdate_url, n_max = 1, progress = FALSE),
    error = function(e) character()
  )
  
  if (length(procdate) == 1 && grepl("^[0-9]{8}$", procdate)) {
    return(procdate)
  }
  
  stop("Could not determine latest nClimDiv processing date from procdate.txt")
}

ensure_latest_pcp_file <- function(cache_dir = CACHE_DIR, base_url = CLIMDIV_BASE_URL, prefix = PCP_PREFIX) {
  ensure_cache_dir(cache_dir)
  
  procdate <- get_latest_procdate(base_url)
  local_path <- file.path(cache_dir, paste0(prefix, "-", procdate))
  remote_url <- paste0(base_url, "/", prefix, "-", procdate)
  
  if (!file.exists(local_path)) {
    utils::download.file(remote_url, destfile = local_path, mode = "wb", quiet = TRUE)
  }
  
  normalizePath(local_path, winslash = "/", mustWork = TRUE)
}

ensure_latest_roni_file <- function(cache_dir = CACHE_DIR, force_refresh = TRUE) {
  ensure_cache_dir(cache_dir)
  
  local_path <- file.path(cache_dir, RONI_CACHE_FILE)
  
  if (force_refresh || !file.exists(local_path)) {
    utils::download.file(RONI_URL, destfile = local_path, mode = "wb", quiet = TRUE)
  }
  
  normalizePath(local_path, winslash = "/", mustWork = TRUE)
}

normalize_enso_class <- function(x) {
  dplyr::case_when(
    x <= -0.5 ~ "La Nina",
    x >=  0.5 ~ "El Nino",
    TRUE ~ "Neutral"
  )
}

# -----------------------------
# Shapefile / lookup
# -----------------------------
read_division_sf <- function(data_dir = DATA_DIR, shapefile_name = SHAPEFILE_NAME) {
  shp_path <- file.path(data_dir, shapefile_name)
  
  if (!file.exists(shp_path)) {
    stop(glue("Shapefile not found: {shp_path}"))
  }
  
  sf::read_sf(shp_path) |>
    dplyr::transmute(
      state_code = as.integer(STATE_CODE),
      state_name = as.character(STATE),
      state_abbr = as.character(ST_ABBRV),
      state_fips = as.character(STATE_FIPS),
      division = as.integer(CD_NEW),
      division_chr = stringr::str_pad(as.integer(CD_NEW), 2, pad = "0"),
      division_name = as.character(NAME),
      climdiv_code = as.integer(CLIMDIV),
      fips_cd = as.character(FIPS_CD),
      geometry = geometry
    )
}

build_division_lookup <- function(div_sf) {
  div_sf |>
    sf::st_drop_geometry() |>
    dplyr::distinct() |>
    dplyr::arrange(state_name, division)
}

# -----------------------------
# Fast nClimDiv precipitation parser
# Wide format: one row per state/division/year with m1:m12
# -----------------------------
parse_climdiv_precip_wide <- function(file_path) {
  dt <- data.table::fread(
    file_path,
    header = FALSE,
    sep = " ",
    strip.white = TRUE,
    fill = TRUE,
    data.table = FALSE,
    showProgress = FALSE,
    colClasses = c("character", rep("numeric", 12))
  )
  
  if (ncol(dt) < 13) {
    stop("Unexpected precipitation file format: fewer than 13 columns.")
  }
  
  dt <- dt[, 1:13]
  names(dt) <- c("id_block", paste0("m", 1:12))
  
  # Ensure fixed width with leading zeros preserved
  dt$id_block <- stringr::str_pad(dt$id_block, width = 10, side = "left", pad = "0")
  
  tibble::as_tibble(dt) |>
    dplyr::transmute(
      state_code = as.integer(substr(id_block, 1, 2)),
      division   = as.integer(substr(id_block, 3, 4)),
      element_code = substr(id_block, 5, 6),
      year       = as.integer(substr(id_block, 7, 10)),
      dplyr::across(dplyr::starts_with("m"), as.numeric)
    ) |>
    dplyr::filter(element_code == "01") |>
    dplyr::select(-element_code)
}

# -----------------------------
# Seasonal precipitation from wide monthly table
# Returns long seasonal summary: state_code, division, year, seasonal_precip_in
# -----------------------------
make_seasonal_precip <- function(monthly_df, season) {
  if (is.null(monthly_df) || nrow(monthly_df) == 0) {
    return(tibble::tibble())
  }
  
  months <- SEASON_MONTHS[[season]]
  if (is.null(months)) {
    stop(glue("Unknown season: {season}"))
  }
  
  # Seasons fully within a calendar year
  if (!identical(months, c(12, 1, 2)) && !identical(months, c(11, 12, 1))) {
    month_cols <- paste0("m", months)
    
    return(
      monthly_df |>
        dplyr::transmute(
          state_code,
          division,
          year,
          seasonal_precip_in = rowSums(dplyr::across(dplyr::all_of(month_cols)), na.rm = TRUE)
        )
    )
  }
  
  # DJF: Dec of previous year + Jan/Feb of current year
  if (identical(months, c(12, 1, 2))) {
    prev_dec <- monthly_df |>
      dplyr::transmute(
        state_code,
        division,
        year = year + 1L,
        dec_prev = m12
      )
    
    curr_jf <- monthly_df |>
      dplyr::transmute(
        state_code,
        division,
        year,
        jan = m1,
        feb = m2
      )
    
    return(
      curr_jf |>
        dplyr::inner_join(prev_dec, by = c("state_code", "division", "year")) |>
        dplyr::transmute(
          state_code,
          division,
          year,
          seasonal_precip_in = dec_prev + jan + feb
        )
    )
  }
  
  # NDJ: Nov/Dec of current year + Jan of next year, assign to next year
  if (identical(months, c(11, 12, 1))) {
    curr_nd <- monthly_df |>
      dplyr::transmute(
        state_code,
        division,
        year = year + 1L,
        nov = m11,
        dec = m12
      )
    
    next_j <- monthly_df |>
      dplyr::transmute(
        state_code,
        division,
        year,
        jan = m1
      )
    
    return(
      next_j |>
        dplyr::inner_join(curr_nd, by = c("state_code", "division", "year")) |>
        dplyr::transmute(
          state_code,
          division,
          year,
          seasonal_precip_in = nov + dec + jan
        )
    )
  }
  
  tibble::tibble()
}

# -----------------------------
# CPC RONI parser
# -----------------------------
get_roni_data <- function(file_path) {
  txt <- readr::read_lines(file_path, progress = FALSE)
  
  if (length(txt) == 0) {
    stop("RONI feed is empty.")
  }
  
  dat <- utils::read.table(
    text = txt,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  names(dat) <- tolower(names(dat))
  
  required <- c("seas", "yr", "anom")
  if (!all(required %in% names(dat))) {
    stop("RONI feed does not contain expected columns: SEAS YR ANOM")
  }
  
  tibble::as_tibble(dat) |>
    dplyr::transmute(
      season = toupper(trimws(seas)),
      year = as.integer(yr),
      roni = as.numeric(anom)
    ) |>
    dplyr::filter(
      !is.na(year),
      !is.na(roni),
      season %in% names(SEASON_MONTHS)
    ) |>
    dplyr::mutate(
      enso_class = normalize_enso_class(roni)
    )
}

# -----------------------------
# Assembly helpers
# -----------------------------
make_plot_dataset <- function(monthly_precip, lookup, roni_df, state_code, division, season) {
  seasonal_precip <- make_seasonal_precip(monthly_precip, season)
  
  if (nrow(seasonal_precip) == 0 || is.null(roni_df) || nrow(roni_df) == 0) {
    return(tibble::tibble())
  }
  
  seasonal_precip |>
    dplyr::filter(state_code == !!state_code, division == !!division) |>
    dplyr::left_join(lookup, by = c("state_code", "division")) |>
    dplyr::inner_join(
      roni_df |>
        dplyr::filter(season == !!season) |>
        dplyr::select(year, season, roni, enso_class),
      by = "year"
    ) |>
    dplyr::arrange(year)
}

# -----------------------------
# UI
# -----------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .control-label { font-weight: 600; }
      .small-note { color: #555; font-size: 0.92em; }
      .metric-card { background: #f7f7f7; border: 1px solid #ddd; border-radius: 8px; padding: 12px 14px; margin-bottom: 12px; }
      .metric-value { font-size: 1.4em; font-weight: 700; }
      .leaflet-container { background: #ffffff; }
    "))
  ),
  titlePanel("ENSO–Precipitation Explorer for U.S. Climate Divisions"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("state", "State", choices = NULL),
      uiOutput("division_ui"),
      selectInput("season", "Season", choices = names(SEASON_MONTHS), selected = "DJF"),
      sliderInput("year_range", "Years", min = 1950, max = 2026, value = c(1950, 2026), sep = ""),
      checkboxInput("show_labels", "Show year labels", FALSE),
      checkboxInput("show_fit", "Show regression line", FALSE),
      checkboxInput("show_ref_lines", "Show reference lines", TRUE),
      br(),
      downloadButton("download_plot", "Download PNG"),
      downloadButton("download_data", "Download CSV"),
      br(), br(),
      div(
        class = "small-note",
        "The app uses NOAA nClimDiv monthly divisional precipitation and CPC Relative Oceanic Niño Index (RONI).",
        tags$br()
      ),
      tags$hr(),
      div(
        style = "font-size: 0.9em; color: #444;",
        tags$b("Contact:"), tags$br(),
        "Mike Crimmins", tags$br(),
        "University of Arizona", tags$br(),
        tags$a(href = "mailto:crimmins@arizona.edu", "crimmins@arizona.edu")
      )
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(4, div(class = "metric-card", textOutput("metric_n_label"), div(class = "metric-value", textOutput("metric_n")))),
        column(4, div(class = "metric-card", textOutput("metric_corr_label"), div(class = "metric-value", textOutput("metric_corr")))),
        column(4, div(class = "metric-card", textOutput("metric_mean_label"), div(class = "metric-value", textOutput("metric_mean"))))
      ),
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          "Scatterplot",
          plotlyOutput("scatter_plot", height = 760) |> shinycssloaders::withSpinner()
        ),
        tabPanel(
          "Data",
          br(),
          tableOutput("plot_data_tbl")
        ),
        tabPanel(
          "Map",
          br(),
          leafletOutput("division_map", height = 700) |> shinycssloaders::withSpinner()
        ),
        tabPanel(
          "About",
          br(),
          verbatimTextOutput("about_txt")
        )
      )
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  div_sf <- reactiveVal(NULL)
  div_lookup <- reactiveVal(NULL)
  monthly_precip <- reactiveVal(NULL)
  roni_data <- reactiveVal(NULL)
  precip_file_used <- reactiveVal(NA_character_)
  
  observeEvent(TRUE, {
    message("Reading climate division shapefile...")
    t0 <- Sys.time()
    
    sf_obj <- read_division_sf()
    lookup <- build_division_lookup(sf_obj)
    
    div_sf(sf_obj)
    div_lookup(lookup)
    
    state_choices <- lookup |>
      dplyr::distinct(state_code, state_name, state_abbr) |>
      dplyr::arrange(state_name) |>
      dplyr::mutate(label = paste0(state_name, " (", state_abbr, ")"))
    
    updateSelectInput(
      session,
      "state",
      choices = stats::setNames(state_choices$state_code, state_choices$label),
      selected = 2
    )
    
    message(glue("Shapefile loaded in {elapsed_secs(t0)} seconds."))
  }, once = TRUE)
  
  observe({
    req(div_lookup())
    
    state_code <- suppressWarnings(as.integer(input$state %||% 2))
    
    div_choices <- div_lookup() |>
      dplyr::filter(state_code == !!state_code) |>
      dplyr::arrange(division) |>
      dplyr::mutate(label = paste0("Division ", division, ": ", division_name))
    
    req(nrow(div_choices) > 0)
    
    sel <- if (state_code == 2 && any(div_choices$division == 1, na.rm = TRUE)) {
      1
    } else {
      div_choices$division[1]
    }
    
    updateSelectInput(
      session,
      "division",
      choices = stats::setNames(div_choices$division, div_choices$label),
      selected = sel
    )
  })
  
  observeEvent(TRUE, {
    message("Preparing cache directory...")
    ensure_cache_dir()
    
    tryCatch({
      message("Checking latest nClimDiv precipitation file...")
      t0 <- Sys.time()
      
      local_or_downloaded <- ensure_latest_pcp_file()
      precip_file_used(basename(local_or_downloaded))
      
      message(glue(
        "Precip file ready: {basename(local_or_downloaded)} ({elapsed_secs(t0)} seconds)."
      ))
      message("Parsing divisional precipitation records...")
      
      t1 <- Sys.time()
      monthly_precip(parse_climdiv_precip_wide(local_or_downloaded))
      message(glue("Precipitation parsed in {elapsed_secs(t1)} seconds."))
    }, error = function(e) {
      message("Precip load failed: ", e$message)
      monthly_precip(NULL)
    })
    
    tryCatch({
      message("Refreshing CPC RONI file...")
      t2 <- Sys.time()
      
      roni_file <- ensure_latest_roni_file(force_refresh = TRUE)
      message(glue(
        "RONI file ready: {basename(roni_file)} ({elapsed_secs(t2)} seconds)."
      ))
      message("Parsing RONI values...")
      
      t3 <- Sys.time()
      roni_data(get_roni_data(roni_file))
      message(glue("RONI parsed in {elapsed_secs(t3)} seconds."))
    }, error = function(e) {
      message("RONI load failed: ", e$message)
      roni_data(NULL)
    })
    
    if (!is.null(monthly_precip()) && !is.null(roni_data())) {
      message("Startup data loading complete.")
    }
  }, once = TRUE)
  
  observe({
    req(roni_data())
    req(nrow(roni_data()) > 0)
    
    yrs <- sort(unique(roni_data()$year))
    
    updateSliderInput(
      session,
      "year_range",
      min = min(yrs),
      max = max(yrs),
      value = c(min(yrs), max(yrs))
    )
  })
  
  selected_info <- reactive({
    req(div_lookup(), input$state, input$division)
    
    div_lookup() |>
      dplyr::filter(
        state_code == as.integer(input$state),
        division == as.integer(input$division)
      ) |>
      dplyr::slice(1)
  })
  
  plot_df <- reactive({
    req(monthly_precip(), div_lookup(), roni_data(), input$state, input$division, input$season)
    req(nrow(monthly_precip()) > 0)
    req(nrow(roni_data()) > 0)
    
    df <- make_plot_dataset(
      monthly_precip = monthly_precip(),
      lookup = div_lookup(),
      roni_df = roni_data(),
      state_code = as.integer(input$state),
      division = as.integer(input$division),
      season = input$season
    )
    
    req(nrow(df) > 0)
    
    df |>
      dplyr::filter(year >= input$year_range[1], year <= input$year_range[2]) |>
      dplyr::mutate(
        enso_class = factor(as.character(enso_class), levels = ENSO_LEVELS)
      )
  })
  
  seasonal_map_summary <- reactive({
    req(monthly_precip(), roni_data(), div_sf(), input$season)
    req(nrow(monthly_precip()) > 0)
    req(nrow(roni_data()) > 0)
    
    season_df <- make_seasonal_precip(monthly_precip(), input$season)
    req(nrow(season_df) > 0)
    
    season_df |>
      dplyr::inner_join(
        roni_data() |>
          dplyr::filter(season == input$season) |>
          dplyr::select(year, roni, enso_class),
        by = "year"
      ) |>
      dplyr::group_by(state_code, division) |>
      dplyr::summarise(
        n_years = dplyr::n(),
        corr = suppressWarnings(cor(roni, seasonal_precip_in, use = "complete.obs")),
        mean_precip = mean(seasonal_precip_in, na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  output$division_ui <- renderUI({
    req(div_lookup(), input$state)
    
    div_choices <- div_lookup() |>
      dplyr::filter(state_code == as.integer(input$state)) |>
      dplyr::arrange(division)
    
    req(nrow(div_choices) > 0)
    
    selectInput(
      "division",
      "Climate division",
      choices = stats::setNames(
        div_choices$division,
        paste0("Division ", div_choices$division, ": ", div_choices$division_name)
      )
    )
  })
  
  output$scatter_plot <- plotly::renderPlotly({
    df <- plot_df()
    req(nrow(df) > 1)
    
    info <- selected_info()
    req(nrow(info) == 1)
    
    corr_val <- suppressWarnings(cor(df$roni, df$seasonal_precip_in, use = "complete.obs"))
    avg_all <- mean(df$seasonal_precip_in, na.rm = TRUE)
    avg_el  <- mean(df$seasonal_precip_in[df$enso_class == "El Nino"], na.rm = TRUE)
    avg_la  <- mean(df$seasonal_precip_in[df$enso_class == "La Nina"], na.rm = TRUE)
    
    df <- df |>
      dplyr::mutate(
        hover_txt = paste0(
          "Year: ", year,
          "<br>RONI: ", sprintf("%.2f", roni),
          "<br>Precipitation: ", sprintf("%.2f in.", seasonal_precip_in),
          "<br>ENSO: ", as.character(enso_class)
        )
      )
    
    p <- ggplot(df, aes(
      x = roni,
      y = seasonal_precip_in,
      color = enso_class,
      text = hover_txt
    )) +
      geom_point(size = 3.6, alpha = 0.95) +
      scale_color_manual(values = ENSO_COLORS, drop = FALSE) +
      labs(
        title = glue("{info$state_name}, Climate Division {info$division}: ENSO vs. Seasonal Precipitation"),
        subtitle = info$division_name,
        x = "Seasonal Average RONI",
        y = "Seasonal Total Precipitation (in.)",
        color = NULL
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = c(0.88, 0.12),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "black"),
        panel.grid.minor = element_blank()
      ) +
      annotate(
        "label",
        x = -Inf,
        y = Inf,
        hjust = -0.05,
        vjust = 1.1,
        size = 4.2,
        label.size = 0.4,
        label.padding = unit(0.25, "lines"),
        fontface = "bold",
        label = glue(
          "Relative Oceanic Nino Index\n",
          "{input$season}\n",
          "{min(df$year)}-{max(df$year)}\n",
          "corr = {round(corr_val, 3)}"
        )
      )
    
    if (isTRUE(input$show_ref_lines)) {
      p <- p +
        geom_hline(yintercept = avg_all, color = "black", linewidth = 0.6) +
        geom_hline(yintercept = avg_el, color = ENSO_COLORS[["El Nino"]], linewidth = 0.6) +
        geom_hline(yintercept = avg_la, color = ENSO_COLORS[["La Nina"]], linewidth = 0.6)
    }
    
    # if (isTRUE(input$show_labels)) {
    #   p <- p + ggrepel::geom_text_repel(
    #     aes(label = year),
    #     size = 4.1,
    #     show.legend = FALSE,
    #     max.overlaps = Inf
    #   )
    # }
    
    # gp <- plotly::ggplotly(p, tooltip = "text") |>
    #   plotly::layout(
    #     legend = list(
    #       x = 0.82,
    #       y = 0.15,
    #       bgcolor = "rgba(255,255,255,0.85)"
    #     )
    #   )
    
    gp <- plotly::ggplotly(p, tooltip = "text") |>
      plotly::layout(
        legend = list(
          x = 1.02,
          y = 0.5,
          xanchor = "left",
          yanchor = "middle",
          bgcolor = "rgba(255,255,255,0.85)",
          bordercolor = "black",
          borderwidth = 1
        ),
        margin = list(r = 140)
      )
    
    if (isTRUE(input$show_fit)) {
      fit <- stats::lm(seasonal_precip_in ~ roni, data = df)
      
      x_seq <- seq(min(df$roni, na.rm = TRUE), max(df$roni, na.rm = TRUE), length.out = 200)
      fit_df <- data.frame(
        roni = x_seq,
        seasonal_precip_in = stats::predict(fit, newdata = data.frame(roni = x_seq))
      )
      
      gp <- gp |>
        plotly::add_lines(
          data = fit_df,
          x = ~roni,
          y = ~seasonal_precip_in,
          inherit = FALSE,
          line = list(color = "gray30", width = 2),
          name = "Linear fit",
          hoverinfo = "skip",
          showlegend = FALSE
        )
    }
    
    if (isTRUE(input$show_labels)) {
      gp <- gp |>
        plotly::add_text(
          data = df,
          x = ~roni,
          y = ~seasonal_precip_in,
          text = ~as.character(year),
          inherit = FALSE,
          textposition = "top center",
          showlegend = FALSE,
          hoverinfo = "skip"
        )
    }
    
    gp
  })
  
  output$plot_data_tbl <- renderTable({
    plot_df() |>
      dplyr::transmute(
        year,
        season,
        roni = round(roni, 2),
        enso_class = as.character(enso_class),
        seasonal_precip_in = round(seasonal_precip_in, 2),
        state = state_name,
        division,
        division_name
      )
  }, striped = TRUE, bordered = TRUE, hover = TRUE, digits = 2)
  
  output$division_map <- renderLeaflet({
    req(input$main_tabs == "Map")
    req(div_sf())
    req(nrow(div_sf()) > 0)
    
    message("Rendering map...")
    t0 <- Sys.time()
    
    map_sf <- div_sf() |>
      dplyr::left_join(seasonal_map_summary(), by = c("state_code", "division"))
    
    pal <- leaflet::colorNumeric(
      palette = "RdBu",
      domain = map_sf$corr,
      reverse = TRUE,
      na.color = "#d9d9d9"
    )
    
    out <- leaflet(map_sf) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addPolygons(
        layerId = ~paste(state_code, division, sep = "_"),
        weight = ~ifelse(
          state_code == as.integer(input$state) & division == as.integer(input$division),
          2.5, 0.8
        ),
        color = ~ifelse(
          state_code == as.integer(input$state) & division == as.integer(input$division),
          "black", "#666666"
        ),
        fillColor = ~pal(corr),
        fillOpacity = 0.75,
        smoothFactor = 0.2,
        label = ~paste0(
          state_abbr, " Division ", division, ": ", division_name,
          "\nCorrelation: ", ifelse(is.na(corr), "NA", sprintf("%.3f", corr)),
          "\nMean precip: ", ifelse(is.na(mean_precip), "NA", sprintf("%.2f in.", mean_precip))
        )
      ) |>
      addLegend(
        position = "bottomright",
        pal = pal,
        values = ~corr,
        title = glue("{input$season} corr"),
        opacity = 0.9
      )
    
    message(glue("Map rendered in {elapsed_secs(t0)} seconds."))
    out
  })
  
  observeEvent(input$division_map_shape_click, {
    click <- input$division_map_shape_click
    req(click$id)
    
    parts <- strsplit(click$id, "_")[[1]]
    
    if (length(parts) == 2) {
      updateSelectInput(session, "state", selected = as.integer(parts[1]))
      updateSelectInput(session, "division", selected = as.integer(parts[2]))
    }
  })
  
  output$metric_n_label <- renderText("Years")
  output$metric_corr_label <- renderText("Correlation")
  output$metric_mean_label <- renderText("Mean precip")
  
  output$metric_n <- renderText({
    df <- tryCatch(plot_df(), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return("0")
    as.character(nrow(df))
  })
  
  output$metric_corr <- renderText({
    df <- tryCatch(plot_df(), error = function(e) NULL)
    if (is.null(df) || nrow(df) < 2) return("NA")
    sprintf("%.3f", suppressWarnings(cor(df$roni, df$seasonal_precip_in, use = "complete.obs")))
  })
  
  output$metric_mean <- renderText({
    df <- tryCatch(plot_df(), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return("NA")
    sprintf("%.2f in", mean(df$seasonal_precip_in, na.rm = TRUE))
  })
  
  output$about_txt <- renderText({
    paste(
      "OVERVIEW
This tool explores the relationship between ENSO (El Niño–Southern Oscillation) and seasonal precipitation across U.S. climate divisions using NOAA datasets.

Each point in the scatterplot represents a single year, showing how seasonal precipitation relates to ENSO conditions during that same period.


DATA SOURCES
- ENSO index: CPC Relative Oceanic Niño Index (RONI)
- Precipitation: NOAA nClimDiv monthly climate division dataset

RONI is a 3-month running mean of Niño 3.4 SST anomalies:
  El Niño  >=  0.5
  La Niña  <= -0.5
  Neutral  between -0.5 and 0.5


HOW TO USE

1. Select a state and climate division
   - Climate divisions represent regions with similar climate behavior

2. Choose a season
   - Seasons are 3-month periods (e.g., DJF = Dec–Feb)
   - Winter seasons span calendar years (e.g., Dec 1999 + Jan–Feb 2000)

3. Adjust year range if needed

4. Toggle display options:
   - Year labels
   - Regression line
   - Reference lines


INTERPRETING THE SCATTERPLOT

- X-axis: ENSO strength (RONI)
- Y-axis: seasonal total precipitation (inches)

Each point:
  One year of data

Colors:
  Blue  = La Niña
  Green = Neutral
  Red   = El Niño

Key features:

• Correlation (corr)
  Positive  = wetter during El Niño
  Negative  = wetter during La Niña
  Near zero = weak relationship

• Regression line
  Shows overall linear trend

• Horizontal lines
  Black = overall mean precipitation
  Red   = El Niño mean
  Blue  = La Niña mean


MAP TAB

The map shows correlation between ENSO and precipitation for each climate division:
  Red shades  = positive correlation
  Blue shades = negative correlation


TECHNICAL NOTES

- Static shapefile assets are read from data/
- Downloaded NOAA files are cached in cache/
- Precipitation is parsed in wide form for speed
- Seasonal totals are computed directly from monthly columns
- Only complete 3-month seasons are included


LIMITATIONS

- Linear correlation does not capture nonlinear ENSO effects
- Relationships may vary over time
- Small sample sizes for short periods can reduce reliability
",
      sep = "\n"
    )
  })
  
  output$download_data <- downloadHandler(
    filename = function() {
      info <- selected_info()
      glue("enso_precip_{info$state_abbr}_div{info$division}_{input$season}.csv")
    },
    content = function(file) {
      readr::write_csv(plot_df(), file)
    }
  )
  
  output$download_plot <- downloadHandler(
    filename = function() {
      info <- selected_info()
      glue("enso_precip_{info$state_abbr}_div{info$division}_{input$season}.png")
    },
    content = function(file) {
      df <- plot_df()
      info <- selected_info()
      
      corr_val <- suppressWarnings(cor(df$roni, df$seasonal_precip_in, use = "complete.obs"))
      avg_all <- mean(df$seasonal_precip_in, na.rm = TRUE)
      avg_el  <- mean(df$seasonal_precip_in[df$enso_class == "El Nino"], na.rm = TRUE)
      avg_la  <- mean(df$seasonal_precip_in[df$enso_class == "La Nina"], na.rm = TRUE)
      
      p <- ggplot(df, aes(x = roni, y = seasonal_precip_in, color = enso_class)) +
        geom_point(size = 3.6, alpha = 0.95) +
        scale_color_manual(values = ENSO_COLORS, drop = FALSE) +
        labs(
          title = glue("{info$state_name}, Climate Division {info$division}: ENSO vs. Seasonal Precipitation"),
          subtitle = info$division_name,
          x = "Seasonal Average RONI",
          y = "Seasonal Total Precipitation (in.)",
          color = NULL
        ) +
        theme_bw(base_size = 13) +
        theme(
          plot.title = element_text(face = "bold"),
          legend.position = c(0.88, 0.12),
          legend.background = element_rect(fill = alpha("white", 0.9), color = "black"),
          panel.grid.minor = element_blank()
        ) +
        annotate(
          "label",
          x = -Inf,
          y = Inf,
          hjust = -0.05,
          vjust = 1.1,
          size = 4.2,
          label.size = 0.4,
          label.padding = unit(0.25, "lines"),
          fontface = "bold",
          label = glue(
            "Relative Oceanic Nino Index\n",
            "{input$season}\n",
            "{min(df$year)}-{max(df$year)}\n",
            "corr = {round(corr_val, 3)}"
          )
        )
      
      if (isTRUE(input$show_ref_lines)) {
        p <- p +
          geom_hline(yintercept = avg_all, color = "black", linewidth = 0.6) +
          geom_hline(yintercept = avg_el, color = ENSO_COLORS[["El Nino"]], linewidth = 0.6) +
          geom_hline(yintercept = avg_la, color = ENSO_COLORS[["La Nina"]], linewidth = 0.6)
      }
      
      if (isTRUE(input$show_fit)) {
        p <- p + geom_smooth(method = "lm", se = FALSE, color = "gray30", linewidth = 0.8)
      }
      
      if (isTRUE(input$show_labels)) {
        p <- p + ggrepel::geom_text_repel(
          aes(label = year),
          size = 4.1,
          show.legend = FALSE,
          max.overlaps = Inf
        )
      }
      
      ggplot2::ggsave(file, plot = p, width = 12, height = 9, dpi = 300)
    }
  )
}

shinyApp(ui, server)