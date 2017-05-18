#' Title
#'
#' @param InputList 
#'
#' @return
#' @export
#'
#' @examples
CompareAcrossData <- function(InputList) {
  
  library(shiny)
  
  # Define UI 
  ui <- fluidPage(
    
    # Application title
    titlePanel("Explore gene expression over the principal curve"),
    
    fluidPage(
      fluidRow(
        wellPanel(
          selectizeInput(
            'GeneList', 'Select gene',
            choices = c('', sort(unique(unlist(lapply(lapply(InputList, "[[", "Expression"), rownames))))),
            selected = '',
            multiple = FALSE)
        ),
        plotOutput("distPlot", height = "1500px")
      )
    )
    
  )
  
  # Define server logic
  server <- function(input, output, session) {
    
    # Functio to ploject on the circle ----------------------------------------
    
    ProjectOnCircle <- function(Points, Edges, Nodes, Categories, Title, ExpValues, PCACenter = TRUE) {
      
      if(PCACenter){
        ScaledNodes <- scale(Nodes, center = PCACenter, scale = FALSE)
        Centers <- attr(ScaledNodes, "scaled:center")
        Nodes <- ScaledNodes
      }
      
      PCAPrGraph <- prcomp(Nodes, retx = TRUE, center = FALSE, scale. = FALSE)
      VarExp <- PCAPrGraph$sdev[1:2]^2/sum(PCAPrGraph$sdev^2)
      
      if(is.null(Categories)){
        Categories <- rep("NoG", nrow(Points))
      }
      
      if(PCACenter){
        Points <- scale(Points, center = Centers, scale = FALSE)
      }
      
      RotatedPoints <- Points %*% PCAPrGraph$rotation[,1:2]
      
      RotatedData <- cbind(RotatedPoints,
                           as.character(Categories),
                           ExpValues[match(rownames(RotatedPoints), names(ExpValues))])
      
      
      colnames(RotatedData) <- c("PC1", "PC2", "Cat", "Exp")
      
      RotatedData.DF <- data.frame(RotatedData)
      RotatedData.DF$PC1 <- as.numeric(as.character(RotatedData.DF$PC1))
      RotatedData.DF$PC2 <- as.numeric(as.character(RotatedData.DF$PC2))
      RotatedData.DF$Exp <- as.numeric(as.character(RotatedData.DF$Exp))
      
      p <- ggplot2::ggplot(data.frame(RotatedData.DF), ggplot2::aes(x = PC1, y = PC2, colour = Exp, shape = Cat)) +
        ggplot2::scale_color_continuous(low = "blue", high = "red")
      
      p <- p + ggplot2::geom_point(data = data.frame(PCAPrGraph$x[,1:2]), mapping = ggplot2::aes(x=PC1, y=PC2), inherit.aes = FALSE) +
        ggplot2::labs(title = Title, x = paste("PC1 -", signif(100*VarExp[1], 4), "%"), y = paste("PC2 -", signif(100*VarExp[2], 4), "%"))
      
      for(j in 1:nrow(Edges)){
        p <- p + ggplot2::geom_path(data = data.frame(PCAPrGraph$x[Edges[j,],1:2]),
                                    mapping = ggplot2::aes(x = PC1, y = PC2), inherit.aes = FALSE)
      }
      
      p <- p + ggplot2::geom_point()
      
      return(p)
      
    }
    
    # Plotting function ----------------------------------------
    
    output$distPlot <- renderPlot({
      
      gName <- input$GeneList
      
      print(paste("Plotting", gName))
      
      if(gName == ''){
        return(NULL)
      }
      
      plotList <- list()
      
      for(i in 1:length(InputList)){
        
        WorkStruct <- InputList[[i]]$PrinCurveStruct$ListProc[[length(InputList[[i]]$PrinCurveStruct$ListProc)]]
        
        ReOrd.Sel <- match(colnames(WorkStruct$CellExp), names(WorkStruct$CellsPT))
        ReOrd.Sel.Mat <- match(colnames(InputList[[i]]$Expression), names(WorkStruct$CellsPT))
        
        if(gName %in% rownames(WorkStruct$NodesExp)){
          
          p1 <- ggplot2::ggplot(data = data.frame(x=WorkStruct$NodesPT,
                                                  y=WorkStruct$NodesExp[gName,]),
                                mapping = ggplot2::aes(x = x, y = y, color="PC")) +
            ggplot2::labs(x = "Pseudotime", y="Gene expression", title = InputList[[i]]$Name) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::geom_rect(data = WorkStruct$RecCoord,
                               mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                               ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
            ggplot2::geom_smooth(data = data.frame(x=WorkStruct$CellsPT[ReOrd.Sel],
                                                   y=WorkStruct$CellExp[gName,]),
                                 mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE) +
            ggplot2::geom_point(data = data.frame(x=WorkStruct$CellsPT[ReOrd.Sel],
                                                  y=WorkStruct$CellExp[gName,]),
                                mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
            ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c(Data = "blue", PC = "black"))
        } else {
          
          if(gName %in% rownames(InputList[[i]]$Expression)){
          
            p1 <- ggplot2::ggplot(data = data.frame(x=WorkStruct$CellsPT[ReOrd.Sel.Mat],
                                                    y=unlist(InputList[[i]]$Expression[gName,])),
                                  mapping = ggplot2::aes(x=x, y=y, color="Data")) +
              ggplot2::labs(x = "Pseudotime", y="Gene expression", title = InputList[[i]]$Name) +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::geom_rect(data = WorkStruct$RecCoord,
                                 mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                                 ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
              ggplot2::geom_smooth() + 
              ggplot2::geom_point(alpha=.5) + ggplot2::scale_color_manual(values = c(Data = "blue", PC = "black"))
            
          } else {
            
            p1 <- ggplot2::ggplot(data = data.frame(x=WorkStruct$CellsPT[ReOrd.Sel.Mat],
                                                    y=unlist(InputList[[i]]$Expression[1,])),
                                  mapping = ggplot2::aes(x=x, y=y, color="Data")) +
              ggplot2::labs(x = "Pseudotime", y="Gene expression", title = InputList[[i]]$Name) +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::geom_rect(data = WorkStruct$RecCoord,
                                 mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                                 ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
              ggplot2::scale_color_manual(values = c(Data = "blue", PC = "black"))
            
          }
          
          
          
        }
        
        
        WorkStruct2 <- InputList[[i]]$PrinCurveStruct$ListInfo[[length(InputList[[i]]$PrinCurveStruct$ListInfo)]]$FinalStruct
        
        if(gName %in% rownames(InputList[[i]]$Expression)){
          p2 <- ProjectOnCircle(
            Points = WorkStruct2$Data,
            Edges = WorkStruct2$IntGrahs[[length(WorkStruct2$IntGrahs)]]$Edges,
            Nodes = WorkStruct2$IntGrahs[[length(WorkStruct2$IntGrahs)]]$Nodes,
            Categories = WorkStruct2$Categories,
            Title = InputList[[i]]$Name,
            ExpValues = unlist(InputList[[i]]$Expression[gName,])
          )
        } else {
          p2 <- ggplot2::ggplot(data = data.frame(x=WorkStruct$NodesPT,
                                                  y=WorkStruct$NodesExp[1,]),
                                mapping = ggplot2::aes(x = x, y = y, color="PC"))
        }
        
        plotList[[length(plotList)+1]] <- p2
        plotList[[length(plotList)+1]] <- p1
        
      }
      
      gridExtra::grid.arrange(grobs = plotList, ncol=2)
      
    })
    
  }
  
  
  # Run the application 
  shinyApp(ui = ui, server = server)
  
}

