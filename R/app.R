#' Visualization bio data through UI
#'
#' bio_shiny function allows user to adjust the visualization of your bio data through the ui
#'
#' @return UI
#' @examples
#'     bio_shiny()
#' @export

bio_shiny <- function() {

setwd("C:/Users/USER/Desktop/r save")

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(shiny)) install.packages("shiny")
if(!require(shinyjs)) install.packages("shinyjs")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggrepel)) install.packages("ggrepel")
if(!require(RColorBrewer)) install.packages("RColorBrewer")

library(ggplot2)
library(shinyjs)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(shiny)

ui <- navbarPage('Test App',id = "inTabset",
                 tabPanel(title = NULL, value = "panel1",
                          fileInput("file_input", "upload pheno type data", buttonLabel = "click"),
                          fileInput("file_input_geno", "upload geno type data", buttonLabel = "click"),
                          fileInput("file_input_adjust", "upload geno type adjust data", buttonLabel = "click"),
                          fileInput("file_input_kendall", "upload kendall data", buttonLabel = "click"),
                          fileInput("file_input_t", "upload t test data", buttonLabel = "click"),
                          ### fileInput("file_input_gwas", "upload gwas data", buttonLabel = "click"),
                          actionButton('jump1to2', 'pheno hist'),
                          actionButton('jump1to3', 'pheno box plot'),
                          actionButton('jump1to4', 'pheno acc dis'),
                          br(),
                          actionButton('jump1to5', 'pheno scatter plot'),
                          actionButton('jump1to6', 'pheno & geno scatter plot'),
                          br(),
                          actionButton('jump1to7', 'pheno & geno box plot'),
                          actionButton('jump1to8', 'pheno & adjust geno scatter plot'),
                          br(),
                          actionButton('jump1to9', 'qq line'),
                          actionButton('jump1to10', 'manhattan plot')),
                 tabPanel(title = NULL, value = "panel2",
                          actionButton('jump2to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              sliderInput("bins",
                                          "Number of bins:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 0.5),
                              selectInput("box_line_color1",
                                          "box line color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("box_size1", "box line size",
                                          choices = seq(0.25, 3.0, 0.25)),
                              selectInput("box_fill_color1",
                                          "box fill color",
                                          choices = c("white", "black", "blue", "red", "yellow", "gray")),
                              selectInput("box_line_type1",
                                          "box line type",
                                          choices = c("solid", "twodash", "longdash", "dotted", "dotdash", "dashed", "blank"))
                            ),

                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot1")
                            )
                          )),
                 tabPanel(title = NULL, value = "panel3",
                          actionButton('jump3to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              sliderInput("transparency2",
                                          "box fill transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("line_length2", "line length",
                                          choices = seq(0.5, 2, 0.5),
                                          selected = 1.5),
                              selectInput("box_line_color2",
                                          "box line color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("box_size2", "box line size",
                                          choices = seq(0.25, 3.0, 0.25)),
                              selectInput("box_fill_color2",
                                          "box fill color",
                                          choices = c("white", "black", "blue", "red", "yellow", "gray")),
                              selectInput("outlier_color2",
                                          "outlier color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("outlier_shape2",
                                          "outlier shape",
                                          choices = c("rectangle" = 0, "circle" = 1, "triangle" = 2, "cross" = 3, "x" = 4, "diamond" = 5)),
                              selectInput("outlier_size2", "outlier size",
                                          choices = seq(1, 10, 1),
                                          selected = 5),
                              selectInput("outlier_stroke2", "outlier strike",
                                          choices = seq(0.25, 2.0, 0.25)),
                              sliderInput("outlier_alpha2",
                                          "outlier transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0)
                            ),
                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot2")
                            )
                          )
                 ),
                 tabPanel(title = NULL, value = "panel4",
                          actionButton('jump4to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("line_color3",
                                          "line color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("line_size3",
                                          "line size",
                                          choices = seq(0.25, 3, 0.25)),
                              sliderInput("transparency3",
                                          "line transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("line_type3",
                                          "line type",
                                          choices = c("solid", "twodash", "longdash", "dotted", "dotdash", "dashed", "blank")),
                              selectInput("line_end_type3",
                                          "line end type",
                                          choices = c("round", "square"))
                            ),

                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot3")
                            )
                          )),
                 tabPanel(title = NULL, value = "panel5",
                          actionButton('jump5to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("point_color4",
                                          "point color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("point_size4",
                                          "point size",
                                          choices = seq(1, 10, 1),
                                          selected = 5),
                              selectInput("point_shape4",
                                          "point shape",
                                          choices = c("circle" = 1, "rectangle" = 0, "triangle" = 2, "cross" = 3, "x" = 4, "diamond" = 5)),
                              sliderInput("point_alpha4",
                                          "point transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("point_stroke4",
                                          "point stroke",
                                          choices = seq(0.25, 3, 0.25))
                            ),

                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot4")
                            )
                          )),
                 tabPanel(title = NULL, value = "panel6",
                          actionButton('jump6to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("point_color5",
                                          "point color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("point_size5",
                                          "point size",
                                          choices = seq(1, 10, 1),
                                          selected = 5),
                              selectInput("point_shape5",
                                          "point shape",
                                          choices = c("circle" = 1, "rectangle" = 0, "triangle" = 2, "cross" = 3, "x" = 4, "diamond" = 5)),
                              sliderInput("point_alpha5",
                                          "point transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("point_stroke5",
                                          "point stroke",
                                          choices = seq(0.25, 3, 0.25)),
                              selectInput("line_color5",
                                          "line color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("line_size5",
                                          "line size",
                                          choices = seq(1, 5, 1),
                                          selected = 2),
                              sliderInput("line_alpha5",
                                          "line transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("line_type5",
                                          "line type",
                                          choices = c("solid", "twodash", "longdash", "dotted", "dotdash", "dashed", "blank"))
                            ),

                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot5")
                            )
                          )),
                 tabPanel(title = NULL, value = "panel7",
                          actionButton('jump7to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              sliderInput("transparency6",
                                          "box fill transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("line_length6", "line length",
                                          choices = seq(0.5, 2, 0.5),
                                          selected = 1.5),
                              selectInput("box_line_color6",
                                          "box line color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("box_size6", "box line size",
                                          choices = seq(0.25, 3.0, 0.25)),
                              selectInput("box_fill_color6",
                                          "box fill color",
                                          choices = c("white", "black", "blue", "red", "yellow", "gray")),
                              selectInput("outlier_color6",
                                          "outlier color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("outlier_shape6",
                                          "outlier shape",
                                          choices = c("rectangle" = 0, "circle" = 1, "triangle" = 2, "cross" = 3, "x" = 4, "diamond" = 5)),
                              selectInput("outlier_size6", "outlier size",
                                          choices = seq(1, 10, 1),
                                          selected = 5),
                              selectInput("outlier_stroke6", "outlier strike",
                                          choices = seq(0.25, 2.0, 0.25)),
                              sliderInput("outlier_alpha6",
                                          "outlier transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("line_color6",
                                          "line color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("line_size6",
                                          "line size",
                                          choices = seq(1, 5, 1),
                                          selected = 2),
                              sliderInput("line_alpha6",
                                          "line transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("line_type6",
                                          "line type",
                                          choices = c("solid", "twodash", "longdash", "dotted", "dotdash", "dashed", "blank"))

                            ),
                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot6")
                            )
                          )),
                 tabPanel(title = NULL, value = "panel8",
                          actionButton('jump8to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("point_color7",
                                          "point color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("point_size7",
                                          "point size",
                                          choices = seq(1, 10, 1),
                                          selected = 5),
                              selectInput("point_shape7",
                                          "point shape",
                                          choices = c("circle" = 1, "rectangle" = 0, "triangle" = 2, "cross" = 3, "x" = 4, "diamond" = 5)),
                              sliderInput("point_alpha7",
                                          "point transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("point_stroke7",
                                          "point stroke",
                                          choices = seq(0.25, 3, 0.25)),
                              selectInput("line_color7",
                                          "line color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("line_size7",
                                          "line size",
                                          choices = seq(1, 5, 1),
                                          selected = 2),
                              sliderInput("line_alpha7",
                                          "line transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0),
                              selectInput("line_type7",
                                          "line type",
                                          choices = c("solid", "twodash", "longdash", "dotted", "dotdash", "dashed", "blank"))
                            ),

                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot7")
                            )
                          )),
                 tabPanel(title = NULL, value = "panel9",
                          actionButton('jump9to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("point_color8",
                                          "point color",
                                          choices = c("black", "white", "blue", "red", "yellow", "gray")),
                              selectInput("point_size8",
                                          "point size",
                                          choices = seq(1, 5, 1),
                                          selected = 3),
                              sliderInput("point_alpha8",
                                          "point transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 1.0)
                            ),

                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot8")
                            )
                          )),
                 tabPanel(title = NULL, value = "panel10",
                          actionButton('jump10to1', 'main'),
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("point_color9",
                                          "point color",
                                          choices = c("Blues", "Greens", "Greys", "Oranges", "Purples", "Reds")),
                              selectInput("highlight_point_color9",
                                          "highlight point color",
                                          choices = c("red", "black", "white", "blue", "yellow", "gray")),
                              textAreaInput("highlight_SNP9",
                                            "highlight SNP (separate with , )",
                                            width = "200px"),
                              sliderInput("point_alpha9",
                                          "point transparency:",
                                          min = 0.0,
                                          max = 1.0,
                                          value = 0.8)

                            ),

                            # Show a plot of the generated distribution
                            mainPanel(
                              plotOutput("distPlot9")
                            )
                          )),
                 useShinyjs(),
                 tags$head(tags$style(HTML('.navbar-nav a {cursor: default}')))
)

options(shiny.maxRequestSize = 1000*1024^2)

server <- function(input, output, session) {

  shinyjs::hide(selector = '.navbar-nav a')
  shinyjs::disable(selector = '.navbar-nav a')

  observeEvent(input$jump1to2, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel2")
  })
  observeEvent(input$jump1to3, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel3")
  })
  observeEvent(input$jump1to4, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel4")
  })
  observeEvent(input$jump1to5, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel5")
  })
  observeEvent(input$jump1to6, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel6")
  })
  observeEvent(input$jump1to7, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel7")
  })
  observeEvent(input$jump1to8, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel8")
  })
  observeEvent(input$jump1to9, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel9")
  })
  observeEvent(input$jump1to10, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel10")
  })

  observeEvent(input$jump2to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })

  observeEvent(input$jump3to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })

  observeEvent(input$jump4to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })

  observeEvent(input$jump5to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })

  observeEvent(input$jump6to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })

  observeEvent(input$jump7to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })

  observeEvent(input$jump8to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })

  observeEvent(input$jump9to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })

  observeEvent(input$jump10to1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  })


  output$distPlot1 <- renderPlot({
    # generate bins based on input$bins from ui.R
    x    <- read.table(input$file_input$datapath , header = F)
    bins <- input$bins
    box_line_color <- input$box_line_color1
    box_size <- as.numeric(input$box_size1)
    box_fill_color <- input$box_fill_color1
    box_line_type <- input$box_line_type1

    ggplot(x, aes(V1)) + geom_histogram(binwidth = bins, aes(y = ..density.., fill = ..count..), color = box_line_color, size = box_size, fill = box_fill_color, linetype = box_line_type) + geom_density(colour = "red") + scale_fill_gradient(low = "white", high = "white") + xlab("value") + ggtitle("pheno type hist") + theme(legend.position = "none" ,plot.title = element_text(hjust = 0.5, size = 25), axis.title = element_text(size = 20), axis.text = element_text(size = 15))
  })

  output$distPlot2 <- renderPlot({
    # generate bins based on input$bins from ui.R
    x    <- read.table(input$file_input$datapath , header = F)
    transparency <- input$transparency2
    line_length <- as.numeric(input$line_length2)
    box_line_color <- input$box_line_color2
    box_size <- as.numeric(input$box_size2)
    box_fill_color <- input$box_fill_color2
    outlier_color <- input$outlier_color2
    outlier_shape <- as.numeric(input$outlier_shape2)
    outlier_size <- as.numeric(input$outlier_size2)
    outlier_stroke <- as.numeric(input$outlier_stroke2)
    outlier_alpha <- as.numeric(input$outlier_alpha2)

    # draw the histogram with the specified number of bins
    # hist(x, breaks = bins, col = 'darkgray', border = 'white')

    ggplot(x, aes(x = 1, y = V1)) + geom_boxplot(alpha = transparency, coef = line_length, color = box_line_color, size = box_size, fill = box_fill_color, outlier.color = outlier_color, outlier.shape = outlier_shape, outlier.size = outlier_size, outlier.stroke = outlier_stroke, outlier.alpha = outlier_alpha) + ggtitle("pheno type boxplot") + theme(plot.title = element_text(hjust = 0.5, size = 25), axis.title = element_text(size = 20), axis.text.x = element_blank(), axis.text.y = element_text(size = 15), legend.title = element_text(size = 15), axis.ticks.x = element_blank()) + xlab(NULL) + ylab("value")
  })

  output$distPlot3 <- renderPlot({
    pheno    <- read.table(input$file_input$datapath, header = F)
    pheno_x <- sort(pheno[, 1])
    df <- data.frame((x = pheno_x), y = pnorm(pheno_x, mean = 0, sd = 1))
    line_color <- input$line_color3
    line_size <- as.numeric(input$line_size3)
    transparency <- as.numeric(input$transparency3)
    line_type <- input$line_type3
    line_end_type <- input$line_end_type3

    ggplot(df, aes(x = x, y = y)) + geom_line(color = line_color, size = line_size, alpha = transparency, linetype = line_type, lineend = line_end_type) + ggtitle("acc line") + theme(plot.title = element_text(hjust = 0.5, size = 25), axis.title = element_text(size = 20), axis.text = element_text(size = 15), legend.title = element_text(size = 15)) + xlab("x") + ylab("y")
  })

  output$distPlot4 <- renderPlot({
    x    <- read.table(input$file_input$datapath, header = F)
    point_color <- input$point_color4
    point_size <- as.numeric(input$point_size4)
    point_shape <- as.numeric(input$point_shape4)
    point_alpha <- as.numeric(input$point_alpha4)
    point_stroke <- as.numeric(input$point_stroke4)

    ggplot(x, aes(x = 1:length(x$V1), y = x$V1)) + geom_point(color = point_color, size = point_size, shape = point_shape, alpha = point_alpha, stroke = point_stroke) + ggtitle("pheno type scatter plot") + theme(plot.title = element_text(hjust = 0.5, size = 25), axis.title = element_text(size = 20), axis.text = element_text(size = 15), legend.title = element_text(size = 15)) + xlab("index") + ylab("value")

  })

  output$distPlot5 <- renderPlot({
    pheno    <- read.table(input$file_input$datapath, header = F)
    geno <- read.table(input$file_input_geno$datapath, header = F)

    x <- data.frame(geno = as.integer(geno[100, ]), pheno = pheno$V1)

    point_color <- input$point_color5
    point_size <- as.numeric(input$point_size5)
    point_shape <- as.numeric(input$point_shape5)
    point_alpha <- as.numeric(input$point_alpha5)
    point_stroke <- as.numeric(input$point_stroke5)
    line_color <- input$line_color5
    line_size <- as.numeric(input$line_size5)
    line_alpha <- as.numeric(input$line_alpha5)
    line_type <- input$line_type5

    ggplot(x, aes(x = geno, y = pheno)) + geom_point(color = point_color, size = point_size, shape = point_shape, alpha = point_alpha, stroke = point_stroke) + geom_abline(size = line_size, color = line_color, alpha = line_alpha, linetype = line_type) + ggtitle("pheno & geno, scatter plot & line") + theme(plot.title = element_text(hjust = 0.5, size = 25), axis.title = element_text(size = 20), axis.text = element_text(size = 15), legend.title = element_text(size = 15))
  })

  output$distPlot6 <- renderPlot({
    pheno    <- read.table(input$file_input$datapath, header = F)
    geno <- read.table(input$file_input_geno$datapath, header = F)

    x <- data.frame(geno = as.integer(geno[100, ]), pheno = pheno$V1)

    transparency <- input$transparency6
    line_length <- as.numeric(input$line_length6)
    box_line_color <- input$box_line_color6
    box_size <- as.numeric(input$box_size6)
    box_fill_color <- input$box_fill_color6
    outlier_color <- input$outlier_color6
    outlier_shape <- as.numeric(input$outlier_shape6)
    outlier_size <- as.numeric(input$outlier_size6)
    outlier_stroke <- as.numeric(input$outlier_stroke6)
    outlier_alpha <- as.numeric(input$outlier_alpha6)

    line_color <- input$line_color6
    line_size <- as.numeric(input$line_size6)
    line_alpha <- as.numeric(input$line_alpha6)
    line_type <- input$line_type6

    ggplot(x, aes(x = geno, y = pheno, group = as.factor(geno), color = as.factor(geno))) + geom_boxplot(alpha = transparency, coef = line_length, color = box_line_color, size = box_size, fill = box_fill_color, outlier.color = outlier_color, outlier.shape = outlier_shape, outlier.size = outlier_size, outlier.stroke = outlier_stroke, outlier.alpha = outlier_alpha) + geom_abline(size = line_size, color = line_color, alpha = line_alpha, linetype = line_type) + ggtitle("qqline") + theme(plot.title = element_text(hjust = 0.5, size = 25), axis.title = element_text(size = 20), axis.text.x = element_blank(), axis.text.y = element_text(size = 15), legend.title = element_text(size = 15), axis.ticks.x = element_blank()) + xlab(NULL) + ylab("value") + scale_color_manual(values=c('black','red', 'blue'))
  })

  output$distPlot7 <- renderPlot({
    pheno    <- read.table(input$file_input$datapath, header = F)
    adjust <- read.table(input$file_input_adjust$datapath, header = F)

    x <- data.frame(adjust = as.double(adjust[100,]), pheno = pheno$V1)

    point_color <- input$point_color7
    point_size <- as.numeric(input$point_size7)
    point_shape <- as.numeric(input$point_shape7)
    point_alpha <- as.numeric(input$point_alpha7)
    point_stroke <- as.numeric(input$point_stroke7)
    line_color <- input$line_color7
    line_size <- as.numeric(input$line_size7)
    line_alpha <- as.numeric(input$line_alpha7)
    line_type <- input$line_type7

    ggplot(x, aes(x = adjust, y = pheno)) + geom_point(color = point_color, size = point_size, shape = point_shape, alpha = point_alpha, stroke = point_stroke) + geom_abline(size = line_size, color = line_color, alpha = line_alpha, linetype = line_type) + ggtitle("pheno & geno, scatter plot & line") + theme(plot.title = element_text(hjust = 0.5, size = 25), axis.title = element_text(size = 20), axis.text = element_text(size = 15), legend.title = element_text(size = 15))

  })

  output$distPlot8 <- renderPlot({
    kendall <- read.table(input$file_input_kendall$datapath, header = F)
    ttest <- read.table(input$file_input_t$datapath, header = F)

    q_kendall <- quantile(kendall$V3, seq(0, 1, 0.01))
    q_ttest <- quantile(ttest$V3, seq(0, 1, 0.01))

    kendall_ttest_data <- data.frame(kendall = q_kendall, ttest = q_ttest)

    point_color <- input$point_color8
    point_size <- as.numeric(input$point_size8)
    point_alpha <- as.numeric(input$point_alpha8)

    ggplot(kendall_ttest_data, aes(kendall, ttest)) + geom_point(color = point_color, size = point_size, shape = 16, alpha = point_alpha) + ggtitle("qq plot") + theme(plot.title = element_text(hjust = 0.5, size = 25), axis.title = element_text(size = 20), axis.text = element_text(size = 15), legend.title = element_text(size = 15))
  })

  output$distPlot9 <- renderPlot({
    df <- read.table(input$file_input_kendall$datapath, header = F)
    df <- data.frame(CHR = df$V1, SNP = df$V2, P = df$V3, BP = df$V4)

    point_color <- input$point_color9
    highlight_point_color <- input$highlight_point_color9
    highlight_SNP <- strsplit(input$highlight_SNP9, ",")
    point_alpha <- as.numeric(input$point_alpha9)

    mypalette <- brewer.pal(7, point_color)[3:7] # chr color palette
    mysnps <- highlight_SNP[[1]]
    # mysnps <- c("rs1000","rs2000","rs3000") # snps to highlight
    sig = 5e-8 # significant threshold line
    sugg = 1e-6 # suggestive threshold line

    gg.manhattan <- function(df, threshold, hlight, col, ylims, title, axis_x_size, anno_size, anno_force){
      # format df
      df.tmp <- df %>%

        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len=max(BP)) %>%

        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%

        # Add this info to the initial dataset
        left_join(df, ., by=c("CHR"="CHR")) %>%

        # Add a cumulative position of each SNP
        arrange(CHR, BP) %>%
        mutate( BPcum=BP+tot) %>%

        # Add highlight and annotation information
        mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
        mutate( is_annotate=ifelse(P < threshold, "yes", "no"))

      # get chromosome center positions for x-axis
      axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

      ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
        # Show all points
        geom_point(aes(color=as.factor(CHR)), alpha=point_alpha, size=2) +
        scale_color_manual(values = rep(col, 22 )) +

        # custom X axis:
        scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis

        # add plot and axis titles
        ggtitle(paste0(title)) +
        labs(x = "Chromosome") +

        # add genome-wide sig and sugg lines
        geom_hline(yintercept = -log10(sig)) +
        geom_hline(yintercept = -log10(sugg), linetype="dashed") +

        # Add highlighted points
        geom_point(data=subset(df.tmp, is_highlight=="yes"), color=highlight_point_color, size=2) +

        # Add label using ggrepel to avoid overlapping
        geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=anno_size, force=anno_force) +

        # Custom the theme:
        theme_bw(base_size = 22) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = axis_x_size)
        )
    }

    gg.manhattan(df, threshold=1e-6, hlight=mysnps, col=mypalette, ylims=c(0,10), title="Manhattan Plot", axis_x_size = 12, anno_size = 5, anno_force = 15) + theme(axis.text.y = element_text(size = 15, color = "black"), axis.text.x = element_text(color = "black"))

  })

}

shinyApp(ui=ui,server=server)


}

# bio_shiny()
