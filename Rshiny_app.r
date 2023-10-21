library(shiny)
library(edgeR)

	## data preprocessing
	# sample_nam <- c("xxx", "xxx", ...) 
	# cond1 <- c("xxx", "xxx", ...)

	# Meta <- cbind(paste0(cond1, "_", sample_nam), cond1)
	# rownames(Meta) <- sample_nam
	# colnames(Meta) <- c("cond0", "cond1")
	# saveRDS(Meta, "Meta.Rds")

	# Counts.matrix <- read.table("counts.txt", skip = 0, header = TRUE)

	# Groups <- factor(Meta[, "cond1"])

	# y_1 <- DGEList(counts = Counts.matrix[, -c(1,2)], group= Groups, lib.size= colSums(Counts.matrix[, -c(1,2)]), norm.factors = rep(1,ncol(Counts.matrix[, -c(1,2)])), samples = NULL, genes = Counts.matrix[,2], remove.zeros = FALSE) 

	# keep <- keep <- filterByExpr(y_1)

	# y_1 <- y_1[keep, , keep.lib.sizes=FALSE]

	# y_1$samples$lib.size <- colSums(y_1$counts)

	# y_1 <- calcNormFactors(y_1)

	# y_1 <- estimateDisp(y_1, robust=TRUE)

	# saveRDS(y_1, "y_1.Rds")


y_1 <- readRDS("y_1.Rds")

Meta <- readRDS("Meta.Rds")

genenames <- rownames(y_1$counts)

Group_options <- list(factor(Meta[, "cond1"], levels=c("Lock","Unedit")), factor(Meta[, "cond0"]))


yourGraph <- function(yourGene, ylimMin = NULL, ylimMax = NULL, Group_options_num = 1){

	yourCounts <- cpm(y_1)[match(yourGene, rownames(cpm(y_1))), , drop=FALSE]

	Groups <- Group_options[[Group_options_num]]
	  
	yourMeans <- apply(yourCounts, 1, function(x){aggregate(x ~ Groups, data=yourCounts, mean, drop=FALSE)})
	
	yourSD <- apply(yourCounts, 1, function(x){aggregate(x ~ Groups, data=yourCounts, sd, drop=FALSE)})

	if(length(ylimMax) > 0){ylimMaxChoice <- ylimMax} else {
		ylimMaxChoice <- max(yourCounts[match(yourGene, rownames(yourCounts)),]) + 0.2*(max(yourCounts[match(yourGene, rownames(yourCounts)),]))}
	
	if(length(ylimMin) > 0){ylimMinChoice <- ylimMin} else {
		ylimMinChoice <- min(yourCounts[match(yourGene, rownames(yourCounts)),]) - 0.2*(max(yourCounts[match(yourGene, rownames(yourCounts)),]))}
	
	stripchart(yourCounts[match(yourGene, rownames(yourCounts)),] ~ Groups,  vert=T, method= "jitter", jitter=0.05, ylim=c(ylimMinChoice, ylimMaxChoice), xlim=c(0.5,length(levels(Groups))+0.5), cex=1, pch=16, axes=TRUE, col=pals::cols25(25)[1:(length(levels(Groups)))], ylab= paste(yourGene, " (cpm)"), cex.axis = 1, las =2, main=paste(yourGene),add=FALSE)

	boxplot(yourCounts[match(yourGene, rownames(yourCounts)),] ~ Groups, col=rgb(red=0.1, green=0.1, blue=0.1, alpha = 0.1), axes=FALSE, add=TRUE)

	for (i in 1:length(levels(Groups))){
		text(rep(i+0.2,table(Groups)[i]),yourCounts[match(yourGene, rownames(yourCounts)),][Groups == levels(Groups)[i]], rownames(Meta)[Groups == levels(Groups)[i]], cex=0.9, col="black")
		}

	mtext(paste0("mean     ", paste0(unlist(yourMeans[[1]][1]), collapse="    ")), side=1, line=9, cex=1, at=c(length(levels(Groups))/2 +0.5))
	
	mtext(paste0("          ", as.numeric(format(round(unlist(yourMeans[[1]][2]), 2), nsmall = 2)), collapse="    "), side=1, line=10, cex=1, at=c(length(levels(Groups))/2 +0.5))
	
	mtext(paste0("sd      ", paste0(unlist(yourSD[[1]][1]), collapse="  ")), side=1, line=11, cex=1, at=c(length(levels(Groups))/2+0.5))
	
	mtext(paste0("        ", as.numeric(format(round(unlist(yourSD[[1]][2]), 2), nsmall=2)), collapse="    "), side=1, line=12, cex=1, at=c(length(levels(Groups))/2 +0.5))

}


Group_options_for_test <- list("AZIN1 locked vs AZIN1 uneditable" = c("Unedit", "Lock"))


ui <- fluidPage(
  titlePanel(h1("HEK293 homozygous AZIN1 locked vs uneditable",align = "center")),
  sidebarLayout(
    sidebarPanel(
      
      selectInput(inputId = "genes",
                  label = "Gene Symbol (case insensitive)",
                  choices = genenames,
                  selected = genenames[3558]),
      
      
      radioButtons(inputId = "whichMeta", label = h3("split by"),
                   choices = list("clone" = 1,
                                  "individual sample" = 2
                   ), 
                   selected = 1),
      
      hr(),
      
      checkboxInput(inputId = "changeRange", label = strong("change y axis"), value = FALSE),
      
      
      conditionalPanel(condition = "input.changeRange == true",
                       sliderInput(inputId = "y_adjust",
                                   label = "y axis adjustment:",
                                   min = 0, max = 3000, value = c(0, 500), step = 2)
      ),
      
      hr(),	
      
	  
	  
	  checkboxInput(inputId = "my_test", label = strong("exact test"), value = FALSE),
	  
      conditionalPanel(condition = "input.my_test == true",
                       radioButtons(inputId = "whichCond", label = h4("Numerator vs Denominator"),
                   choices = list("AZIN1 locked vs AZIN1 uneditable" = 1
                   ), 
                   selected = 1)
      ),
	  



	conditionalPanel(condition = "input.my_test == true",
                       numericInput(inputId = "top_n", label = h4("show top N"), value = 10),
					   fluidRow(column(4, verbatimTextOutput("top_n")))
      ),
	  	  
      conditionalPanel(condition = "input.my_test == true",
	  fluidRow(column(12, verbatimTextOutput("my_test_result")))
	  
	  ),	


      width = 3),
    
    mainPanel(
      fluidRow(
        style= "border: 10px dashed white",splitLayout(
          cellWidths = c("100%"),
          plotOutput(outputId = "main_plot",height = '800px',width = "1100px")
        )
      )
    )
  ),
  

  
  hr(),
  fluidRow(
    column(width = 1, align = "center", img(src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/47/Indiana_Hoosiers_logo.svg/150px-Indiana_Hoosiers_logo.svg.png", width='100%')),
    column(width = 11,
           p(
              "XXXXX", a("XXXXX", href="XXXXXXXXXXXXXXX")
           )
    )
  ),
  br()
)



server <- function(input, output, session) {
  
  observe({if(input$changeRange == TRUE){
    output$main_plot <- renderPlot({ yourGraph(input$genes, 
                                               ylimMin = input$y_adjust[1], ylimMax = input$y_adjust[2], 
                                               Group_options_num = as.numeric(input$whichMeta)) })
  } else {
    output$main_plot <- renderPlot({ yourGraph(input$genes, 
                                               ylimMin = NULL, ylimMax = NULL, 
                                               Group_options_num = as.numeric(input$whichMeta))
    })
  }
  })					
  

  output$whichCond <- renderPrint({ input$whichCond })	
  
  observe({if(input$my_test == TRUE){
  output$my_test_result <- renderPrint({ 
				topTags(exactTest(y_1, 
				pair = c(Group_options_for_test[[as.numeric(input$whichCond)]][1],Group_options_for_test[[as.numeric(input$whichCond)]][2])), 
				n=as.numeric(input$top_n)) })
		}
	})

}


shinyApp(ui = ui, server = server)

