library(shiny)

# analytic PDF for x>0 using Bessel I
d_compound_exp <- function(x, lambda, mu, t) {
  xpos <- x > 0
  dens <- rep(0, length(x))
  if (any(xpos)) {
    z <- lambda * t * mu * x[xpos]
    dens[xpos] <- exp(-lambda*t - mu*x[xpos]) *
      sqrt((lambda*t*mu) / x[xpos]) *
      besselI(2*sqrt(z), 1)
  }
  dens
}

# simulate S(t)
r_compound_exp <- function(nsim, lambda, mu, t) {
  N <- rpois(nsim, lambda*t)
  S <- numeric(nsim)
  uniq <- unique(N)
  
  for (k in uniq) {
    idx <- which(N == k)
    if (k == 0) S[idx] <- 0
    else S[idx] <- rgamma(length(idx), shape = k, rate = mu)
  }
  S
}

ui <- fluidPage(
  
  titlePanel("Compound Poisson with Exponential Jumps: S(t)"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("lambda", "Inter-arrival rate λ:", 0.01, 5, 0.5, step = 0.01),
      sliderInput("mu", "Jump rate μ:", 0.01, 5, 1, step = 0.01),
      numericInput("nsim", "Number of simulations:", 5000, min = 100),
      numericInput("seed", "Random seed:", 2025, min = 0),
      checkboxInput("show_density", "Overlay analytic density (x>0)", TRUE)
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("t=10", plotOutput("p10"), verbatimTextOutput("i10")),
        tabPanel("t=100", plotOutput("p100"), verbatimTextOutput("i100")),
        tabPanel("t=1000", plotOutput("p1000"), verbatimTextOutput("i1000")),
        tabPanel("t=10000", plotOutput("p10000"), verbatimTextOutput("i10000"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  simulateS <- function(t) {
    set.seed(input$seed)
    r_compound_exp(input$nsim, input$lambda, input$mu, t)
  }
  
  renderPlotFor <- function(S, t) {
    lam <- input$lambda
    mu  <- input$mu
    
    # histogram safe limits
    h <- hist(S, breaks = 50, plot = FALSE)
    ymax <- max(h$density, na.rm = TRUE)
    if (!is.finite(ymax)) ymax <- 1  # fallback
    
    hist(S, breaks = 50, freq = FALSE,
         col = "lightblue",
         main = paste("S(t) at t =", t),
         xlab = "S(t)",
         ylim = c(0, ymax))
    
    if (input$show_density) {
      xgrid <- seq(1e-6, max(S) + 1, length.out = 600)
      lines(xgrid, d_compound_exp(xgrid, lam, mu, t), col = "red", lwd = 2)
    }
  }
  
  renderInfoFor <- function(S, t) {
    lam <- input$lambda
    mu  <- input$mu
    
    paste0(
      "t = ", t, "\n",
      "Simulated Mean = ", mean(S), "\n",
      "Simulated Var  = ", var(S), "\n",
      "Theoretical Mean = ", lam*t/mu, "\n",
      "Theoretical Var  = ", 2*lam*t/mu^2, "\n",
      "P(S=0) = exp(-λt) = ", exp(-lam*t), "\n"
    )
  }
  
  output$p10 <- renderPlot({ renderPlotFor(simulateS(10), 10) })
  output$i10 <- renderText({ renderInfoFor(simulateS(10), 10) })
  
  output$p100 <- renderPlot({ renderPlotFor(simulateS(100), 100) })
  output$i100 <- renderText({ renderInfoFor(simulateS(100), 100) })
  
  output$p1000 <- renderPlot({ renderPlotFor(simulateS(1000), 1000) })
  output$i1000 <- renderText({ renderInfoFor(simulateS(1000), 1000) })
  
  output$p10000 <- renderPlot({ renderPlotFor(simulateS(10000), 10000) })
  output$i10000 <- renderText({ renderInfoFor(simulateS(10000), 10000) })
}

shinyApp(ui,server)
