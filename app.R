###############################################################
# FULL R SHINY APP — Compound Poisson Process with Exp Jumps
# S(t) = sum_{i=1}^{N(t)} X_i
# N(t) ~ Poisson(λ t)
# X_i ~ Exp(μ)
###############################################################

# ------------------------
# Load required packages
# ------------------------
install.packages("shiny")
library(shiny)
library(ggplot2)

# ------------------------------------------------------------
# THEORETICAL DENSITY OF S(t)
# S(t) has mass e^{-λt} at 0. For y>0:
# f(y) = e^{-λt - μy} * sqrt(λtμ / y) * I1(2 * sqrt(λtμy))
# ------------------------------------------------------------

dCompoundExp <- function(y, lambda, mu, t) {
  dens <- numeric(length(y))
  idx_pos <- which(y > 0)
  
  if (length(idx_pos) > 0) {
    z <- 2 * sqrt(lambda * t * mu * y[idx_pos])
    b <- besselI(z, 1)
    dens[idx_pos] <- exp(-lambda * t - mu * y[idx_pos]) *
      sqrt(lambda * t * mu / y[idx_pos]) * b
  }
  
  return(dens)
}

# ------------------------------------------------------------
# SIMULATION FUNCTION
# ------------------------------------------------------------
rCompoundExp <- function(nsim, lambda, mu, t) {
  N <- rpois(nsim, lambda * t)
  S <- numeric(nsim)
  
  idx <- which(N > 0)
  for (i in idx) {
    S[i] <- sum(rexp(N[i], rate = mu))
  }
  return(S)
}

# ------------------------------------------------------------
# THEORETICAL MOMENTS
# ------------------------------------------------------------
meanCompound <- function(lambda, mu, t) lambda * t / mu
varCompound  <- function(lambda, mu, t) 2 * lambda * t / (mu^2)


# ------------------------------------------------------------
# USER INTERFACE
# ------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Compound Poisson Process with Exponential Jumps"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("lambda", "Arrival rate (λ):", min = 0.01, max = 5,
                  value = 0.5, step = 0.01),
      sliderInput("mu", "Jump rate (μ):", min = 0.01, max = 5,
                  value = 1, step = 0.01),
      sliderInput("t", "Time t:", min = 1, max = 20000,
                  value = 10, step = 1),
      numericInput("nsim", "Number of simulations:",
                   value = 20000, min = 100, max = 300000),
      checkboxInput("overlay", "Overlay theoretical density", TRUE),
      sliderInput("bins", "Histogram bins:", min = 20, max = 200,
                  value = 80, step = 5),
      actionButton("resim", "Run Simulation")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Single t",
                 plotOutput("histPlot", height = "500px"),
                 verbatimTextOutput("stats")
        ),
        tabPanel("Multiple t (10, 100, 1000, 10000)",
                 plotOutput("manyHist", height = "800px"),
                 helpText("Shows histograms at four fixed time points.")
        )
      )
    )
  )
)

# ------------------------------------------------------------
# SERVER LOGIC
# ------------------------------------------------------------
server <- function(input, output, session) {
  
  # Run simulation reactively
  sim <- eventReactive(input$resim, {
    rCompoundExp(input$nsim, input$lambda, input$mu, input$t)
  }, ignoreNULL = FALSE)
  
  # -----------------------------  
  # PLOT: Histogram for chosen t
  # -----------------------------
  output$histPlot <- renderPlot({
    S <- sim()
    df <- data.frame(S = S)
    
    p <- ggplot(df, aes(x = S)) +
      geom_histogram(aes(y = ..density..),
                     bins = input$bins, fill = "grey80", color = "black") +
      labs(title = paste("Histogram of S(t) at t =", input$t),
           x = "S(t)", y = "Density")
    
    # Overlay theoretical density
    if (input$overlay) {
      xmax <- quantile(S, 0.995)
      xgrid <- seq(0, xmax, length.out = 2000)
      dens <- dCompoundExp(xgrid, input$lambda, input$mu, input$t)
      p <- p + geom_line(data = data.frame(x=xgrid, y=dens),
                         aes(x=x, y=y), color = "blue", linewidth = 1)
      
      # mark P(S=0)
      p <- p + annotate(
        "text", x = xmax * 0.05,
        y = max(dens, na.rm = TRUE) * 0.9,
        label = paste0("P(S=0) = ", round(exp(-input$lambda * input$t), 4)),
        hjust = 0, color = "blue"
      )
    }
    p
  })
  
  # -----------------------------  
  # THEORETICAL STATS
  # -----------------------------
  output$stats <- renderPrint({
    cat("Theoretical Results:\n")
    cat("---------------------\n")
    cat("E[S(t)] =", meanCompound(input$lambda, input$mu, input$t), "\n")
    cat("Var[S(t)] =", varCompound(input$lambda, input$mu, input$t), "\n")
    cat("P(S(t) = 0) =", exp(-input$lambda * input$t), "\n")
  })
  
  # -----------------------------  
  # MULTIPLE t VALUES
  # -----------------------------
  output$manyHist <- renderPlot({
    times <- c(10, 100, 1000, 10000)
    par(mfrow = c(2, 2))
    
    for (tt in times) {
      S <- rCompoundExp(20000, input$lambda, input$mu, tt)
      hist(S, breaks = input$bins, probability = TRUE,
           main = paste("t =", tt), xlab = "S(t)", col = "grey90")
      
      if (input$overlay) {
        xgrid <- seq(0, quantile(S, 0.995), length.out = 2000)
        dens <- dCompoundExp(xgrid, input$lambda, input$mu, tt)
        lines(xgrid, dens, col = "blue", lwd = 2)
        mtext(paste("P(S=0) =", round(exp(-input$lambda*tt), 4)),
              side = 3, line = -2, col = "blue")
      }
    }
  })
}

# ------------------------------------------------------------
# RUN THE APP
# ------------------------------------------------------------
shinyApp(ui, server)


runApp("app.R")
