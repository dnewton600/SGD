#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)
library(ggplot2)

grad <- function(a,b,c,theta) {
  return( c( 2*a*theta[1] + 2*c*theta[2], 2*b*theta[2] + 2*c*theta[1] ) )
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("SGD Simulation for 2D Quadratic Form"),
   titlePanel(h4("f(x,y) = ax^2 + by^2 + 2cxy")),
   titlePanel(h4("f'(x,y) = (2ax + 2cy, 2by + 2cx)")),
   titlePanel(h4("error ~ N(0, var = M0 + M1*||f'(x,y)||^2)")),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectInput("select", "Algorithm", 
                    choices = list("SGD: step size = alpha" = 1, 
                                   "SGD: step size = alpha/k" = 2,
                                   "ADAM" = 3,
                                   "Momentum" = 4), 
                    selected = 1),
        numericInput("a", "a:", 1),
        numericInput("b", "b:", 1),
        numericInput("c", "c:", 0),
        numericInput("M0", "M0:", .1),
        numericInput("M1", "M1:", .1),
        numericInput("alpha", "alpha:", .2),
        actionButton(inputId = "go", label = "Run")
      ), 
      mainPanel(plotOutput("my_hist")))
)

# Define server logic 
server <- function(input, output) {
  
  data <- eventReactive(input$go, {
    
    # get points of the unit circle
    x <- seq(-1,1, by = .01)
    y_upper <- sqrt(1 - x^2)
    y_lower <- -sqrt(1 - x^2)
    
    #plot(x = c(x,-x), y = c(y_lower, y_upper), asp = 1, type = 'l', xlim = c(-2,2), ylim = c(-2,2))
    
    my_mat <- matrix( c(1/input$a, -input$c, -input$c, 1/input$b), nrow = 2)
    
    points <- rbind(c(x,-x), c(y_lower,y_upper))
    
    points_transformed <- my_mat%*%points
    plot(x = points_transformed[1,], y = points_transformed[2,],
         xlab = 'x', ylab = 'y', asp = 1, type = 'l', xlim = c(-3,3), ylim = c(-3,3))
    lines(x = .5*points_transformed[1,], y = .5*points_transformed[2,], type = 'l')
    lines(x = 1.5*points_transformed[1,], y = 1.5*points_transformed[2,], type = 'l')
    lines(x = 2*points_transformed[1,], y = 2*points_transformed[2,], type = 'l')
    
    theta_init <- c(2,1)
    theta <- theta_init
    
    num_iters <- 20
    record <- matrix(0, nrow = num_iters, ncol = 2)
    
    #adam parameters
    beta_1 <- 0.9
    beta_2 <- 0.999
    epsilon_adam <- .000001
    m <- 0
    v <- 0
    
    # momentum parameters
    moment <- .1
    
    for(ii in 1:num_iters) {
      true_grad <- grad(input$a, input$b, input$c, theta)
      sd <- sqrt(input$M0 + input$M1*sum( true_grad^2 ))
      noise_norm <- rnorm(1, sd = sd)
      noise <- rnorm(2)
      noise_unit_vec <- noise/(sqrt(sum(noise^2)))
      
      if(input$select == 1) { 
        theta <- theta - input$alpha*(true_grad + noise_unit_vec*noise_norm)
      }
      
      else if(input$select == 2) {
        theta <- theta - (input$alpha/ii)*(true_grad + noise_unit_vec*noise_norm)
      }
      
      else if(input$select == 3) {
        
        g <- true_grad + noise_unit_vec*noise_norm
        m <- beta_1*m + (1 - beta_1)*g
        v <- beta_2*v + (1 - beta_2)*g^2
        m_hat <- m/(1 - beta_1^ii)
        v_hat <- v/(1 - beta_2^ii)
        theta <- theta - input$alpha*(m_hat/(sqrt(v_hat) + epsilon_adam))
      }
      
      else if(input$select == 4) {
        # theta_prev <- ifelse(ii==1,theta,record[ii-1,])
        # theta <- theta - (input$alpha/ii)*(true_grad + noise_unit_vec*noise_norm) + moment*(theta - theta_prev)
        
        g <- true_grad + noise_unit_vec*noise_norm
        m <- beta_1*m + (1 - beta_1)*g
        m_hat <- m/(1 - beta_1^ii)
        theta <- theta - input$alpha*m_hat
      }
      
      record[ii,] <- theta
    }
    
    record <- rbind(matrix(theta_init, nrow = 1), record)
    
    lines(x = record[,1], y = record[,2], type = 'l', col = 'red')
    
  })
  
  output$my_hist <- renderPlot({
    data()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

