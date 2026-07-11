#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

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
      helpText("Adjust the following conditions to observe the bias due to the scoring strategy and a full latent approach. True value = 0.60."),
      helpText("Adjust the number of items used to model the latent variable:"),
      selectInput("n_items", "Number of Items (Indicators):", 
                  choices = c(3, 6, 12, 18), selected = 3),
      helpText("Adjust the measurement quality of the items (similar to reliability):"),
      selectInput("communality", "Item Communality:", 
                  choices = c("High", "Low"), selected = 'High'),
      helpText("Adjust the sample size:"),
      selectInput("N", "Sample Size:", 
                  choices = c(30, 60, 100, 200, 400, 500, 1000), selected = 30)
      
    ),
    
    mainPanel(
      plotOutput("mcPlot"),
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
    communalities <- data.frame(LOW = c(0.2, 0.3, 0.4, 0.2, 0.3, 0.4), 
                                #WIDE = c(0.2, 0.32, 0.44, 0.56, 0.68, 0.8),
                                HIGH = c(0.6, 0.7, 0.8, 0.6, 0.7, 0.8))
    phi <- matrix(c(
      1.000, 0.6, # Y
      0.6, 1.000 # X
    ), ncol = 2, byrow = TRUE)
    
    n_items <- as.numeric(input$n_items)
    lambda <- matrix(0, nrow = (n_items*2), ncol = 2)
    idx_Y <- 1:n_items
    num_var <- n_items*2
    idx_X <- (n_items+1):num_var
    #lambda2_X <- 1
    # Set parameters based on communality choice
    if (input$communality == "High") {
      lambda2_X<- lambda2_Y <- rep_len(pull(communalities['HIGH']), n_items)
      lambda[idx_Y, 1] <- sqrt(lambda2_Y)
      lambda[idx_X, 2] <- sqrt(lambda2_X)
      Lambda_full <- lambda
      theta <- diag(num_var) - diag(diag(Lambda_full %*% phi %*% t(Lambda_full)))
      # Implied Cov matrix
      sigma0_kappa <- Lambda_full %*% phi %*% t(Lambda_full) + theta
      sigma0_kappa <- (sigma0_kappa + t(sigma0_kappa))/2 
      colnames(sigma0_kappa) <- c(paste0('y', 1:n_items), paste0('x', 1:n_items))
      
    } else {
      lambda2_X<- lambda2_Y <- rep_len(pull(communalities['LOW']), n_items)
      lambda[idx_Y, 1] <- sqrt(lambda2_Y)
      lambda[idx_X, 2] <- sqrt(lambda2_X)
      Lambda_full <- lambda
      theta <- diag(num_var) - diag(diag(Lambda_full %*% phi %*% t(Lambda_full)))
      # Implied Cov matrix
      sigma0_kappa <- Lambda_full %*% phi %*% t(Lambda_full) + theta
      sigma0_kappa <- (sigma0_kappa + t(sigma0_kappa))/2 
      colnames(sigma0_kappa) <- c(paste0('y', 1:n_items), paste0('x', 1:n_items))
    }
    
    # True Causal Structure (SEM)
    dat <- as.data.frame(rmvnorm(N,
                                 mean = rep(0,num_var),
                                 sigma = sigma0_kappa)) 
    colnames(dat) <- c(paste0('y', 1:n_items), paste0('x', 1:n_items))
    
    model <-  c(paste0("etaY =~ NA*", paste0(paste0("y", 1:n_items), collapse = " + ")),
                paste0("etaX =~ NA*", paste0(paste0("x", 1:n_items), collapse = " + ")),
                "etaY  ~~ varetaY*etaY",
                "etaX  ~~ 1*etaX",
                "etaY  ~ a*etaX",
                "varetaY==1-a^2"
    )
    
    model.fit<-sem(model, data =dat )
    model.params <- parameterEstimates(model.fit, standardized = TRUE)
    b1 <- model.params[which(model.params$label==c('a')), 'est']
    
    # Sum Scores
    vsY <- t(rep(1, n_items)) %*% sigma0_kappa[idx_Y,idx_Y] %*% rep(1, n_items)
    dat$sumY <- 1/sqrt(vsY[1,1])*rowSums(dat[idx_Y])
    vsX <- t(rep(1, n_items)) %*% sigma0_kappa[idx_X,idx_X] %*% rep(1, n_items)
    dat$sumX <- 1/sqrt(vsX[1,1])*rowSums(dat[idx_X])
    
    dat$sumY_std <- as.numeric(scale(dat$sumY))
    dat$sumX_std <- as.numeric(scale(dat$sumX))
    model.sum <- lm(sumY_std ~ sumX_std, data = dat)
    
    # FS Regression & Bartlett
    model.fy <-  c(paste0("FY =~ NA*", paste0(paste0("y", 1:n_items), collapse = " + ")),
                   "FY  ~~ 1*FY"
    )
    model.fx <-  c(paste0("FX =~ NA*", paste0(paste0("x", 1:n_items), collapse = " + ")),
                   "FX  ~~ 1*FX"
    )
    model.fy.fit<-cfa(model.fy, data =dat )
    model.fx.fit<-cfa(model.fx, data =dat )
    dat[, c("FrY")] <- lavPredict(model.fy.fit, method="regression")
    dat[, c("FbY")] <- lavPredict(model.fy.fit, method='Bartlett')
    dat[, c("FrX")] <- lavPredict(model.fx.fit, method="regression")
    dat[, c("FbX")] <- lavPredict(model.fx.fit, method='Bartlett')
    #model.rr <- lm(FrY ~ FrX, data = dat)
    #model.br <- lm(FbY ~ FbX, data = dat) 
    #model.br.r <- lm(FbY ~ FrX, data = dat) 
    
    model.reg <-  c(paste0("etaY =~ NA*", paste0(paste0("y", 1:n_items), collapse = " + ")),
                    "etaY  ~~ varetaY*etaY",
                    "etaY  ~ a*FrX",
                    "varetaY==1-a^2"
    )
    model.reg.fit<-sem(model.reg, data =dat )
    model.reg.params <- parameterEstimates(model.reg.fit, standardized = TRUE)
    b1.reg <- model.params[which(model.reg.params$label==c('a')), 'est']   
    
    model.bar <-  c(paste0("etaY =~ NA*", paste0(paste0("y", 1:n_items), collapse = " + ")),
                    "etaY  ~~ varetaY*etaY",
                    "etaY  ~ a*FbX",
                    "varetaY==1-a^2"
    )
    model.bar.fit<-sem(model.bar, data =dat )
    model.bar.params <- parameterEstimates(model.bar.fit, standardized = TRUE)
    b1.bar <- model.params[which(model.bar.params$label==c('a')), 'est'] 
    
    model.sum <-  c(paste0("etaY =~ NA*", paste0(paste0("y", 1:n_items), collapse = " + ")),
                    "etaY  ~~ varetaY*etaY",
                    "etaY  ~ a*sumX",
                    "varetaY==1-a^2"
    )
    model.sum.fit<-sem(model.sum, data =dat )
    model.sum.params <- parameterEstimates(model.sum.fit, standardized = TRUE)
    b1.sum <- model.params[which(model.sum.params$label==c('a')), 'est'] 
    
    # Estimates table
    df_sem <- data.frame(intercept = 0,  slope = b1, Model = "SEM")
    df_sum <- data.frame(intercept = 0,  slope = b1.sum, Model = "SS")
    df_reg <- data.frame(intercept = 0,  slope = b1.reg, Model = "RS")
    df_bar <- data.frame(intercept = 0,  slope = b1.bar, Model = "BS")
    #df_br <- data.frame(intercept = 0,  slope = b1, Model = "BR-S")
    df_true<- data.frame(intercept = 0,  slope = .6, Model = "True DGP")
    
    return(list(lines_df = rbind(df_sem, df_sum, df_reg, df_bar, df_true), raw_data = dat))
  })
  
  output$mcPlot <- renderPlot({
    res <- sim_data()
    lines_df <- res$lines_df
    raw_data <- res$raw_data
    
    ggplot() +
      # 1. Plot raw data points to establish scales (using x and the calculated sumY)
      geom_point(data = raw_data, aes(x = sumX_std, y = sumY_std), alpha = 0.4, color = "grey") +
      
      # 2. Draw lines using slope and intercept parameters from your simulation dataframe
      geom_abline(data = lines_df, aes(intercept = intercept, slope = slope, color = Model, linetype = Model), linewidth = .8, alpha = 0.5) +
      
      # 3. Aesthetics
      scale_color_manual(values = c("SEM" = "green", "SS" = "pink", "RS" = "#377EB8", "BS" = "#377EB8", "True DGP" = "orange")) +
      scale_linetype_manual(values = c("SEM" = "solid", "SS" = "solid", "RS" = "dotted", "BS" = "dashed", "True DGP" = "solid")) +
      labs(title = "Error free estimation by modeling lantent variables",
           x = "Predictor Variable (X)",
           y = "Outcome Variable (Y)") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$mcSummary <- renderPrint({
    # Add summary readout logic here if desired
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
