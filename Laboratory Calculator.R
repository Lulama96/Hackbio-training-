library(shiny)
library(ggplot2)
library(plotly)

ui <- fluidPage(
 
  tags$style(HTML("
        .btn-custom {
            background-color: #007bff; /* Blue color */
            color: white;
            border: none;
            padding: 10px 20px;
            font-size: 16px;
            border-radius: 4px;
            cursor: pointer;
        }
        .btn-custom:hover {
            background-color: #0056b3; /* Darker blue color on hover */
        }body {
                background-color: #F3E5F5;
            }
            .shiny-output-error {
                color: red;
            }
            .title-box {
                background-color: #4A148C;
                padding: 20px;
                text-align: center;
                border-radius: 10px;
                color: white;
                margin-bottom: 20px;
            }
            .title-box h1 {
                margin: 0;
                font-size: 32px;
            }
    ")),
  
  # Title in Header Box
  div(class = "title-box",
      h1("Laboratory Reagent Prepation Calculators")), 
  sidebarLayout(
    sidebarPanel(
      # Tabs for different calculators
      tabsetPanel(
        id = 'tabs',
        tabPanel("Serial Dilution",
                 numericInput("init_conc", "Initial Concentration (C1)", value = NULL),
                 numericInput("final_conc", "Final Concentration (C2)", value = NULL),
                 numericInput("dilution_factor", "Dilution Factor", value = 2),
                 numericInput("total_volume", "Total Volume", value = NULL),
                 selectInput("volume_unit", "Volume Unit", 
                             choices = c("mL" = "mL", "µL" = "µL", "L" = "L"), selected = "mL"),
                 actionButton("calc_serial", "Calculate Serial Dilution", class = "btn-custom"),
                 actionButton("clear_serial", "Clear", class = "btn-custom")
        ),
        tabPanel("Stock Dilution",
                 numericInput("stock_conc", "Stock Concentration", value = NULL),
                 numericInput("final_conc_stock", "Final Concentration", value = NULL),
                 numericInput("final_volume_stock", "Final Volume", value = NULL),
                 selectInput("volume_unit_stock", "Volume Unit", 
                             choices = c("mL" = "mL", "µL" = "µL", "L" = "L"), selected = "mL"),
                 actionButton("calc_stock", "Calculate Stock Dilution", class = "btn-custom"),
                 actionButton("clear_stock", "Clear", class = "btn-custom")
        ),
        tabPanel("Molarity",
                 numericInput("solutes_molarity", "Molarity (M)", value = NULL),
                 numericInput("solutes_volume", "Volume (L)", value = NULL),
                 numericInput("solutes_mass", "Molecular Weight (g/mol)", value = NULL),
                 actionButton("calc_molarity", "Calculate Mass Required", class = "btn-custom"),
                 actionButton("clear_molarity", "Clear", class = "btn-custom")
        ),
        tabPanel("Buffer Preparation",
                 numericInput("buffer_conc", "Desired Buffer Concentration (M)", value = NULL),
                 numericInput("buffer_volume", "Desired Buffer Volume (L)", value = NULL),
                 numericInput("buffer_stock_conc", "Stock Solution Concentration (M)", value = NULL),
                 actionButton("calc_buffer", "Calculate Volume of Stock Solution", class = "btn-custom"),
                 actionButton("clear_buffer", "Clear", class = "btn-custom")
        )
      )
    ),
    mainPanel(
      # Output area for selected calculator
      textOutput("serial_dilution_output"),
      plotlyOutput("serial_dilution_schematic"),
      textOutput("stock_dilution_output"),
      textOutput("molarity_output"),
      textOutput("buffer_output")
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$calc_serial, {
    C1 <- input$init_conc
    C2 <- input$final_conc
    dilution_factor <- input$dilution_factor
    V2 <- input$total_volume
    
    volume_unit <- input$volume_unit
    conversion_factor <- switch(volume_unit,
                                "mL" = 1,
                                "µL" = 1000,
                                "L" = 0.001)
    V2 <- V2 * conversion_factor
    
    num_steps <- log10(C1 / C2) / log10(dilution_factor)
    volumes <- numeric(num_steps)
    concentrations <- numeric(num_steps)
    
    for (i in 1:num_steps) {
      C_current <- C1 / (dilution_factor ^ i)
      V1 <- (C_current * V2) / C1
      volumes[i] <- V1
      concentrations[i] <- C_current
    }
    volumes <- volumes / conversion_factor
    
    output$serial_dilution_output <- renderText({
      paste("Number of steps: ", round(num_steps),
            "\nVolumes to transfer (in", volume_unit, "): ", paste(round(volumes, 2), collapse = ", "))
    })
    
    output$serial_dilution_schematic <- renderPlotly({
      df <- data.frame(
        step = 1:num_steps,
        concentration = concentrations,
        volume = volumes
      )
      
      rect_width <- 0.8
      rect_height <- 0.5
      
      p <- ggplot(df, aes(x = step, fill = concentration)) +
        geom_rect(aes(xmin = step - rect_width / 2, xmax = step + rect_width / 2,
                      ymin = 0 - rect_height / 2, ymax = 0 + rect_height / 2,
                      text = paste("Step:", step, "<br>Concentration:", round(concentration, 2), 
                                   "<br>Volume:", round(volume, 2), volume_unit)), 
                  color = "black") +
        geom_segment(aes(x = step + rect_width / 2, xend = step + 1 - rect_width / 2,
                         y = 0, yend = 0), 
                     arrow = arrow(length = unit(0.3, "cm")), 
                     size = 1, color = "blue", data = df[df$step < num_steps, ]) +
        scale_fill_gradient(low = "#D8BFD8", high = "#800080", name = "Concentration") +
        theme_void() +
        labs(title = "Serial Dilution Schematic", subtitle = "Test Tubes and Transfers") +
        theme(
          plot.title = element_text(hjust = 0.5, size = 18), 
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid = element_blank()
        )
      
      ggplotly(p, tooltip = "text")
    })
  })
  
  observeEvent(input$calc_stock, {
    stock_conc <- input$stock_conc
    final_conc <- input$final_conc_stock
    final_volume <- input$final_volume_stock
    
    volume_unit_stock <- input$volume_unit_stock
    conversion_factor_stock <- switch(volume_unit_stock,
                                      "mL" = 1,
                                      "µL" = 1000,
                                      "L" = 0.001)
    final_volume <- final_volume * conversion_factor_stock
    
    volume_stock_needed <- (final_conc * final_volume) / stock_conc
    volume_stock_needed <- volume_stock_needed / conversion_factor_stock
    
    output$stock_dilution_output <- renderText({
      paste("Volume of stock solution required (in", volume_unit_stock, "):", round(volume_stock_needed, 2))
    })
  })
  
  observeEvent(input$calc_molarity, {
    molarity <- input$solutes_molarity
    volume <- input$solutes_volume
    molecular_weight <- input$solutes_mass
    
    mass_required <- molarity * volume * molecular_weight
    
    output$molarity_output <- renderText({
      paste("Mass of solute required (in grams):", round(mass_required, 2))
    })
  })
  
  observeEvent(input$calc_buffer, {
    buffer_conc <- input$buffer_conc
    buffer_volume <- input$buffer_volume
    stock_conc <- input$buffer_stock_conc
    
    volume_stock_needed <- (buffer_conc * buffer_volume) / stock_conc
    
    output$buffer_output <- renderText({
      paste("Volume of stock solution required (in L):", round(volume_stock_needed, 2))
    })
  })
  
  observeEvent(input$clear_serial, {
    updateNumericInput(session, "init_conc", value = NULL)
    updateNumericInput(session, "final_conc", value = NULL)
    updateNumericInput(session, "dilution_factor", value = 2)
    updateNumericInput(session, "total_volume", value = NULL)
    updateSelectInput(session, "volume_unit", selected = "mL")
    output$serial_dilution_output <- renderText({ NULL })
    output$serial_dilution_schematic <- renderPlotly({ NULL })
  })
  
  observeEvent(input$clear_stock, {
    updateNumericInput(session, "stock_conc", value = NULL)
    updateNumericInput(session, "final_conc_stock", value = NULL)
    updateNumericInput(session, "final_volume_stock", value = NULL)
    updateSelectInput(session, "volume_unit_stock", selected = "mL")
    output$stock_dilution_output <- renderText({ NULL })
  })
  
  observeEvent(input$clear_molarity, {
    updateNumericInput(session, "solutes_molarity", value = NULL)
    updateNumericInput(session, "solutes_volume", value = NULL)
    updateNumericInput(session, "solutes_mass", value = NULL)
    output$molarity_output <- renderText({ NULL })
  })
  
  observeEvent(input$clear_buffer, {
    updateNumericInput(session, "buffer_conc", value = NULL)
    updateNumericInput(session, "buffer_volume", value = NULL)
    updateNumericInput(session, "buffer_stock_conc", value = NULL)
    output$buffer_output <- renderText({ NULL })
  })
}

shinyApp(ui = ui, server = server)

