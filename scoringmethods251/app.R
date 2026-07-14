#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

conditions_error<-read.csv('https://github.com/Tonatiuhxotl/251-project/blob/main/conditions_error_251app.csv?raw=true')
load(url("https://github.com/Tonatiuhxotl/251-project/blob/main/p_matrix2.rda?raw=true"))
sim251_data<-read.csv('https://github.com/Tonatiuhxotl/251-project/blob/main/sim251_est_values.csv?raw=true')

library(shiny)
library(dplyr)
library(ggplot2)
library(mvtnorm)
library(lavaan)


# Define UI for application that draws a histogram
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h4("Latent Outcome Variable"),
      #selectInput("scenario", "Select Scenario:", 
      #            choices = c("Case 1: Latent Outcome (Y)", 
      #                        "Case 2: Latent Treatment Propensity (X)")),
      hr(),
      h5("Measurement Model Parameters"),
      helpText("Adjust the following conditions to observe the bias due to the scoring strategy"),
      helpText("Adjust the communality level:"),
      selectInput("communality", "Item Communality:", 
                  choices = c("High", "Wide", "Low"), selected = 'High'),
      helpText("Adjust the magnitude of model error (RMSEA) and form (XL or RC):"),
      selectInput("moderror", "Error:", 
                  choices = c("e=0.09 (XL)", "e=0.03 (XL)","e=0", "e=0.03 (RC)", "e=0.09 (RC)"), selected = 'e=0'),
      helpText("Adjust the sample size:"),
      selectInput("N", "Sample Size:", 
                  choices = c(30, 60, 100, 200, 300, 400), selected = 30),
      helpText("Adjust the model specification:"),
      selectInput("specification", "Specification:", 
                  choices = c("Simple", "Correct"), selected = 'Simple')
      
    ),
    
    mainPanel(
      plotOutput("mcPlotA"),
      plotOutput("mcPlotB"),
      plotOutput("mcPlotC"),
      hr()
      #h4("ML Regression Estimates"),
      #verbatimTextOutput("mcSummary")
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  sim_data <- reactive({
    # Fixing seed for stable comparison when parameters are toggled
    set.seed(2026) 
    N <- as.numeric(input$N)
    communality <-as.character(input$communality)
    moderror <- as.character(input$moderror)
    specification <- as.character(input$specification)
    num_var<-18
    n_items<-6
    idx_eta1<-c(1:6)
    idx_eta2<-c(7:12)
    idx_eta3<-c(13:18)
    k<-conditions_error[which(conditions_error$communality == communality & conditions_error$moderror == moderror), 'Id']
    sigma0_kappa <- p_matrix[[k]]
    
    # True Causal Structure (SEM)
    dat <- as.data.frame(rmvnorm(N,
                                 mean = rep(0,num_var),
                                 sigma = )) 
    colnames(dat) <- paste0('y', 1:num_var)
    
    # Sum Scores
    vsEta1 <- t(rep(1, n_items)) %*% sigma0_kappa[idx_eta1,idx_eta1] %*% rep(1, n_items)
    dat$eta1_sum <- 1/sqrt(vsEta1[1,1])*rowSums(dat[idx_eta1])
    vsEta2 <- t(rep(1, n_items)) %*% sigma0_kappa[idx_eta2,idx_eta2] %*% rep(1, n_items)
    dat$eta2_sum <- 1/sqrt(vsEta2[1,1])*rowSums(dat[idx_eta2])
    vsEta3 <- t(rep(1, n_items)) %*% sigma0_kappa[idx_eta3,idx_eta3] %*% rep(1, n_items)
    dat$eta3_sum <- 1/sqrt(vsEta3[1,1])*rowSums(dat[idx_eta3])
    
    dat$eta1_sum_std <- as.numeric(scale(dat$eta1_sum))
    dat$eta2_sum_std <- as.numeric(scale(dat$eta2_sum))
    dat$eta3_sum_std <- as.numeric(scale(dat$eta3_sum))
    
    filter_data <- sim251_data %>% 
      filter(NSIZE == N &
               MODERROR == moderror &
               COMMUNALITY == communality & 
               SPECIFICATION == specification)
    
    a.sem <- filter_data[which(filter_data$PAR=='a' & filter_data$METHOD == 'SEM'), 'est']
    a.ilv <- filter_data[which(filter_data$PAR=='a' & filter_data$METHOD == 'ILV'), 'est']   
    a.sum <- filter_data[which(filter_data$PAR=='a' & filter_data$METHOD == 'SS'), 'est'] 
    a.reg <- filter_data[which(filter_data$PAR=='a' & filter_data$METHOD == 'RS'), 'est'] 
    a.bar <- filter_data[which(filter_data$PAR=='a' & filter_data$METHOD == 'BS'), 'est']
    
    b.sem <- filter_data[which(filter_data$PAR=='b' & filter_data$METHOD == 'SEM'), 'est']
    b.ilv <- filter_data[which(filter_data$PAR=='b' & filter_data$METHOD == 'ILV'), 'est']   
    b.sum <- filter_data[which(filter_data$PAR=='b' & filter_data$METHOD == 'SS'), 'est'] 
    b.reg <- filter_data[which(filter_data$PAR=='b' & filter_data$METHOD == 'RS'), 'est'] 
    b.bar <- filter_data[which(filter_data$PAR=='b' & filter_data$METHOD == 'BS'), 'est']
    
    c.sem <- filter_data[which(filter_data$PAR=='c' & filter_data$METHOD == 'SEM'), 'est']
    c.ilv <- filter_data[which(filter_data$PAR=='c' & filter_data$METHOD == 'ILV'), 'est']   
    c.sum <- filter_data[which(filter_data$PAR=='c' & filter_data$METHOD == 'SS'), 'est'] 
    c.reg <- filter_data[which(filter_data$PAR=='c' & filter_data$METHOD == 'RS'), 'est'] 
    c.bar <- filter_data[which(filter_data$PAR=='c' & filter_data$METHOD == 'BS'), 'est']
    
    # Estimates table
    df_sem <- data.frame(a_path = a.sem, b_path = b.sem,  c_path = c.sem, Model = "SEM")
    df_ilv <- data.frame(a_path = a.ilv, b_path = b.ilv,  c_path = c.ilv, Model = "ILV")
    df_sum <- data.frame(a_path = a.sum, b_path = b.sum,  c_path = c.sum, Model = "SS")
    df_reg <- data.frame(a_path = a.reg, b_path = b.reg,  c_path = c.reg, Model = "RS")
    df_bar <- data.frame(a_path = a.bar, b_path = b.bar,  c_path = c.bar, Model = "BS")
    df_true<- data.frame(a_path = 0.3, b_path = 0.7,  c_path = 0, Model = "True")
    
    return(list(lines_df = rbind(df_sem, df_ilv, df_sum, df_reg, df_bar, df_true), raw_data = dat))
  })
  
  output$mcPlotA <- renderPlot({
    res <- sim_data()
    lines_df <- res$lines_df
    raw_data <- res$raw_data
    
    ggplot() +
      # 1. Plot raw data points to establish scales (using x and the calculated sumY)
      geom_point(data = raw_data, aes(x = eta1_sum_std, y = eta2_sum_std), alpha = 0.4, color = "grey") +
      
      # 2. Draw lines using slope and intercept parameters from your simulation dataframe
      geom_abline(data = lines_df, aes(intercept = 0, slope = a_path, color = Model, linetype = Model), linewidth = .8, alpha = 0.5) +
      
      # 3. Aesthetics
      scale_color_manual(values = c("SEM" = "darkgreen", "ILV" = "darkgreen","SS" = "pink", "RS" = "#377EB8", "BS" = "#377EB8", "True" = "orange")) +
      scale_linetype_manual(values = c("SEM" = "dotted", "ILV" = "dashed", "SS" = "solid", "RS" = "dotted", "BS" = "dashed", "True" = "solid")) +
      labs(title = "A path (b_21)",
           x = "Predictor (Eta1)",
           y = "Mediator (Eta2)") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  
  output$mcPlotB <- renderPlot({
    res <- sim_data()
    lines_df <- res$lines_df
    raw_data <- res$raw_data
    
    ggplot() +
      # 1. Plot raw data points to establish scales (using x and the calculated sumY)
      geom_point(data = raw_data, aes(x = eta2_sum_std, y = eta3_sum_std), alpha = 0.4, color = "grey") +
      
      # 2. Draw lines using slope and intercept parameters from your simulation dataframe
      geom_abline(data = lines_df, aes(intercept = 0, slope = b_path, color = Model, linetype = Model), linewidth = .8, alpha = 0.5) +
      
      # 3. Aesthetics
      scale_color_manual(values = c("SEM" = "darkgreen", "ILV" = "darkgreen","SS" = "pink", "RS" = "#377EB8", "BS" = "#377EB8", "True" = "orange")) +
      scale_linetype_manual(values = c("SEM" = "dotted", "ILV" = "dashed", "SS" = "solid", "RS" = "dotted", "BS" = "dashed", "True" = "solid")) +
      labs(title = "B path (b_32)",
           x = "Mediator (Eta2)",
           y = "Outcome (Eta3)") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$mcPlotC <- renderPlot({
    res <- sim_data()
    lines_df <- res$lines_df
    raw_data <- res$raw_data
    
    ggplot() +
      # 1. Plot raw data points to establish scales (using x and the calculated sumY)
      geom_point(data = raw_data, aes(x = eta1_sum_std, y = eta3_sum_std), alpha = 0.4, color = "grey") +
      
      # 2. Draw lines using slope and intercept parameters from your simulation dataframe
      geom_abline(data = lines_df, aes(intercept = 0, slope = c_path, color = Model, linetype = Model), linewidth = .8, alpha = 0.5) +
      
      # 3. Aesthetics
      scale_color_manual(values = c("SEM" = "darkgreen", "ILV" = "darkgreen","SS" = "pink", "RS" = "#377EB8", "BS" = "#377EB8", "True" = "orange")) +
      scale_linetype_manual(values = c("SEM" = "dotted", "ILV" = "dashed", "SS" = "solid", "RS" = "dotted", "BS" = "dashed", "True" = "solid")) +
      labs(title = "C path (b_31)",
           x = "Predictor (Eta1)",
           y = "Outcome (Eta3)") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$mcSummary <- renderPrint({
    # Add summary readout logic here if desired
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
