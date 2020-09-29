# -----package load and library-----
library(shiny)
library(seqinr)
library(plyr)
library(gdata)
library(ggplot2)

# -----Functions-----

# find the matching point
markRessite <- function(fasta, res_enz, cutting_site) {
  matching_matrix <- gregexpr(res_enz, fasta, ignore.case = TRUE)
  matching_matrix <- data.frame(matching_matrix)
  matched <- as.vector(matching_matrix)[ , 1]
  # output cutting result
  
  under_line <- paste(replicate(matched[1]-1, " "), collapse = "", sep = "")
  lines <- paste(replicate(nchar(res_enz), "-"), collapse = "" , sep = "")
  under_line <- paste(under_line, lines, collapse = "", sep = "")
  
  if (length(matched) >= 2) {
    for (i in 2:length(matched)) {
      spaces <- matched[i]-(matched[i-1]+nchar(res_enz))
      lines <- paste(replicate(nchar(res_enz), "-"), collapse = "" , sep = "")
      append_elem <- paste(replicate(spaces, " "), collapse = "", sep = "")
      under_line <- paste(under_line, append_elem, lines, collapse = "", sep = "")
    }
  }
  
  for (i in 1:length(matched)) {
    under_line <- paste(substr(under_line, 1, matched[i]+cutting_site-3), ">", "<", 
                        substr(under_line, matched[i]+cutting_site, nchar(under_line)), sep = "")
  }
  output <- paste(fasta, "\n", under_line, sep = "")
  return(output)
}  

# return dataframe with where enzymes cut
posRessite <- function(fasta, res_enz, cutting_site) {
  matching_matrix <- gregexpr(res_enz, fasta, ignore.case = TRUE)
  matching_matrix <- data.frame(matching_matrix)
  colnames(matching_matrix) <- "position"
  if (matching_matrix[1,1] == -1) {
    colnames(matching_matrix) <- "position"
    return(matching_matrix)
  }
  return(matching_matrix + cutting_site - 2)
}

# return a data frame with predicted band length after cut
bandRescut <- function(fasta, pos) {
  bandlength <- NULL
  bandlength[1] <- pos[1,1]
  if (length(pos[, 1]) >= 2) {
    for (i in 2:length(pos[, 1])) {
      bandlength[i] <- pos[i,1] - pos[i-1,1]
    }
  }
  bandlength[length(pos[, 1])+1] <- nchar(fasta) - pos[length(pos[, 1]),1]
  return(data.frame(bandlength))
}

# assign restriction enzyme name with the detecting pattern and cutting site
resMatch <- function(enz_name) {
  seq <- switch(enz_name,
                "ApaI" = "GGGCCC",
                "BamHI" = "GGATCC",
                "BglII" = "AGATCT",
                "EcoRI" = "GAATTC",
                "HindIII" = "AAGCTT",
                "KpnI" = "GGTACC",
                "NcoI" = "CCATGG",
                "NdeI" = "CATATG",
                "NheI" = "GCTAGC",
                "NotI" = "GCGGCCGC",
                "SacI" = "GAGCTC",
                "SalI" = "GTCGAC",
                "SphI" = "GCATGC",
                "XbaI" = "TCTAGA",
                "XhoI" = "CTCGAG"
  )
  cut <- switch(enz_name,
                "ApaI" = 6,
                "BamHI" = 2,
                "BglII" = 2,
                "EcoRI" = 2,
                "HindIII" = 2,
                "KpnI" = 6,
                "NcoI" = 2,
                "NdeI" = 3,
                "NheI" = 2,
                "NotI" = 3,
                "SacI" = 6,
                "SalI" = 2,
                "SphI" = 6,
                "XbaI" = 2,
                "XhoI" = 2
  )
  return(c(seq, cut))
}

res_enz_list = c("ApaI","BamHI","BglII","EcoRI","HindIII",
                 "KpnI","NcoI","NdeI","NheI","NotI",
                 "SacI","SalI","SphI","XbaI","XhoI")
# Define UI for application that offer cutting site/predicted band length
ui <- fluidPage(
   
   # Application title
   titlePanel("Find Restriction Enzyme Cutting Site Related To SNP"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         selectInput("restriction_enzyme","Choose a restriction enzyme:", 
                     choices = res_enz_list),
         fileInput('file1', 'Choose FASTA file',
                   accept=c('fasta','fa')),
         helpText("Load the FASTA file contains all the individuals included in RFLP analysis."),
         # submitButton("Update View"),
         
         helpText("Download cutting site table:"),
         downloadButton('downloadData2', 'Download'),
         
         helpText("Download prediction table:"),
         downloadButton('downloadData', 'Download')

         # helpText("Download prediction map:"),
         # downloadButton('downloadData3', 'Download')
      ),


      
      # show 2 tables and 1 predicted gel map
      mainPanel(
         h4("Cutting Site"),
         tableOutput("cutdata"),
         h4("Predicted Band Length"),
         tableOutput("headdata"),
         h4("Predict Gel Map:"),
         plotOutput("gelmap")
      )
   )
)

# Define server logic required to output tables and graphs
server <- function(input, output) {


  # output gel map
  output$gelmap <- renderPlot({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    
    fasta_list <- read.fasta(inFile$datapath, as.string = TRUE)
    res_enz <- resMatch(input$restriction_enzyme)[1]
    cutting_site <- as.numeric(resMatch(input$restriction_enzyme)[2])
    
    pos_table <- posRessite(unlist(fasta_list[1]), res_enz, cutting_site)
    bandpredict_table <- bandRescut(unlist(fasta_list[1]), pos_table)
    if (length(fasta_list) >=2 ) {
      for (i in 2:length(fasta_list)) {
        fasta_iter <- unlist(fasta_list[i])
        pos <- posRessite(fasta_iter, res_enz, cutting_site)
        pos_table <- cbindX(pos_table, posRessite(fasta_iter, res_enz, cutting_site))
        bandpredict_table <- cbindX(bandpredict_table, bandRescut(fasta_iter, pos))
      }
    }
    
    rownames(bandpredict_table) <- paste(rep("band",length(bandpredict_table[,1])), as.character(seq(1,length(bandpredict_table[,1]))))
    
    individual_names <- attr(fasta_list, which = "name")
    
    seq_bandmap <- data.frame(c(NULL,NULL))
    # convert table into plot drawing
    for (i in seq_along(bandpredict_table[1, ])) {
      table_temp <- cbind(data.frame(rep(individual_names[i],length(na.omit(bandpredict_table[, i])))), na.omit(bandpredict_table[,i]))
      seq_bandmap <- rbind(seq_bandmap, table_temp)
    }
    colnames(seq_bandmap) <- c("individuals", "bandlength")
    ladder <- data.frame(c("M","M","M","M","M","M","M","M","M","M","M","M","M","M"),
                         c(50,100,200,300,400,500,600,700,800,900,1000,1500,2000,5000))
    colnames(ladder) <- c("individuals", "bandlength")
    seq_bandmap <- rbind(ladder, seq_bandmap)
    seq_bandmap <- seq_bandmap[which(seq_bandmap$bandlength != -1),]
    
    bandname <- setNames(object = as.character(c("M", individual_names)), nm = c("M", individual_names))
    # match(seq_bandmap$individuals, unique(seq_bandmap$individuals))
    
    plot_gel <- ggplot(seq_bandmap) +
      # scale_x_discrete(name = "individuals", limits = unique(individuals)) +
      theme(axis.text.x = element_text(colour = "black"), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank()) + 
      geom_rect(data = seq_bandmap, aes(xmin = match(individuals, unique(individuals)) - 0.2, 
                                        xmax = match(individuals, unique(individuals)) + 0.2,
                                        ymin = bandlength - 15,
                                        ymax = bandlength + 15))+
      scale_x_discrete(name = "Individual", limits = bandname) +
      scale_y_continuous(trans='reverse', breaks = c(50,100,200,300,400,500,600,700,800,900,1000,1500,2000,5000),
                         labels = (c("50","100","200","300","400","500","600","700","800","900","1000","1500","2000","5000")))
    plot_gel
  })
  # cutting table output  
   output$cutdata <- renderTable({
     inFile <- input$file1
     if (is.null(inFile))
       return(NULL)
     
     fasta_list <- read.fasta(inFile$datapath, as.string = TRUE)
     res_enz <- resMatch(input$restriction_enzyme)[1]
     cutting_site <- as.numeric(resMatch(input$restriction_enzyme)[2])
     
     pos_table <- posRessite(unlist(fasta_list[1]), res_enz, cutting_site)
     for (i in 2:length(fasta_list)) {
       fasta_iter <- unlist(fasta_list[i])
       pos <- posRessite(fasta_iter, res_enz, cutting_site)
       pos_table <- cbindX(pos_table, posRessite(fasta_iter, res_enz, cutting_site))
     }
     rownames(pos_table) <- paste(rep("position",length(pos_table[, 1])), 
                                          as.character(seq(1,length(pos_table[, 1]))))
     pos_table[pos_table == -1] <- NA
     individual_names <- attr(fasta_list, which = "name")
     pos_table <- rbind(individual_names, pos_table)
     t(pos_table)
   })
   # length table output
   output$headdata <- renderTable({
     inFile <- input$file1
     if (is.null(inFile))
       return(NULL)
     
     fasta_list <- read.fasta(inFile$datapath, as.string = TRUE)
     res_enz <- resMatch(input$restriction_enzyme)[1]
     cutting_site <- as.numeric(resMatch(input$restriction_enzyme)[2])
     
     pos_table <- posRessite(unlist(fasta_list[1]), res_enz, cutting_site)
     bandpredict_table <- bandRescut(unlist(fasta_list[1]), pos_table)
     for (i in 2:length(fasta_list)) {
       fasta_iter <- unlist(fasta_list[i])
       pos <- posRessite(fasta_iter, res_enz, cutting_site)
       pos_table <- cbindX(pos_table, posRessite(fasta_iter, res_enz, cutting_site))
       bandpredict_table <- cbindX(bandpredict_table, bandRescut(fasta_iter, pos))
     }
     rownames(bandpredict_table) <- paste(rep("band",length(bandpredict_table[, 1])), 
                                 as.character(seq(1,length(bandpredict_table[, 1]))))
     bandpredict_table[bandpredict_table == -1] <- NA
     individual_names <- attr(fasta_list, which = "name")
     bandpredict_table <- rbind(individual_names, bandpredict_table)
     t(bandpredict_table)
     
   })
# download files  
   output$downloadData2 <- downloadHandler(
     filename = "example1.csv",
     content = function(file) {
       write.csv(pos_table, file)
     })

   output$downloadData <- downloadHandler(
     filename = "example2.csv",
     content = function(file) {
       write.csv(bandpredict_table, file)
     })

   # output$downloadData3 <- downloadHandler(
   #   filename = "example3.png",
   #   content = function(file) {
   #     ggsave(plot_gel, file)
   #   })

   

}

# Run the application 
shinyApp(ui = ui, server = server)

