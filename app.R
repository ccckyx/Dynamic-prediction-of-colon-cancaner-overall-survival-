options(encoding = "UTF-8")
library(animation)
library(ggplot2)
library(scales)
library(ggpubr)
library(funData)
library(MASS)
library(randomForestSRC)
library(pec)
library(survival)
library(fdapace)
library(shiny)
source("./function_all.R")

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Dynamic prediction for overall survival of colorectal cancer patients based on longitudinal CEA, CA19-9, and CA125"),
  
  
  fluidRow(
    fluidRow(
      column(1, offset = 0.5),
      
      column(3,
             numericInput("ID", "ID",
                          value = 311, width = '100%')),
      
      column(3,
             numericInput("NL", "Age",
                          value = 60, min = 18, max = 80, width = '100%')),
      
      column(3,
             selectInput(inputId = "XB",
                         label = "Sex",
                         choices = c("Male", "Female"), width = '100%')),
      
      column(1, offset = 0.5)
    ),
    
    fluidRow(
      column(1, offset = 0.5),
      
      column(3,
             selectInput(inputId = "SSLJ",
                         label = "Surgical approach",
                         choices = c("Laparoscopic resection", "Open resection"), width = '100%')),
      
      column(3,
             selectInput(inputId = "YFBW",
                         label = "Primary site",
                         choices = c("Colon", "Rectum"), width = '100%')),
      
      column(3,
             selectInput(inputId = "FHCD",
                         label = "Tumor differentiation",
                         choices = c("Well", "Moderate", "Poor-undifferentiated", "Unknown"),
                         width = '100%')),
      
      column(1, offset = 0.5)
      
    ),
    
    
    fluidRow(
      column(1, offset = 0.5),
      
      column(3,
             selectInput(inputId = "BLFQ",
                         label = "AJCC 8th ed. Stage",
                         choices = c("I", "II", "III"), width = '100%')),
      
      column(3,
             selectInput(inputId = "LYN",
                         label = "Lymph node yield less than 12",
                         choices = c("Yes", "No"), width = '100%')),
      
      column(3,
             selectInput(inputId = "FZZL",
                         label = "Adjuvant chemotherapy",
                         choices = c("Yes", "No"), width = '100%')),
      column(1, offset = 0.5)
      
    ),
    
    fluidRow(
      column(1, offset = 0.5),
      
      column(3,
             selectInput(inputId = "BLLX",
                         label = "Mucinous (colloid) type",
                         choices = c("Yes", "No"), width = '100%')),
      
      column(3,
             selectInput(inputId = "MGAS",
                         label = "Lymphovascular invasion",
                         choices = c("Yes", "No"), width = '100%')),
      
      column(3,
             selectInput(inputId = "SJSQF",
                         label = "Perineural invasion",
                         choices = c("Yes", "No"), width = '100%')),
      
      column(1, offset = 0.5)
      
    ),
    
    fluidRow(
      column(1, offset = 0.5),
      
      column(3,
             numericInput("JX_CEA", "Preoperative CEA",
                          value = 5, min = 0, max = 1000, width = '100%')),
      
      column(3,
             numericInput("JX_CA199", "Preoperative CA19-9",
                          value = 37, min = 0, max = 1000, width = '100%')),
      
      column(3,
             numericInput("JX_CA125", "Preoperative CA125",
                          value = 35, min = 0, max = 1000, width = '100%')),
      
      column(1, offset = 0.5)
      
    ),
    
    
    fluidRow(
      column(1, offset = 0.5),
      
      column(9,
             fileInput("file1", "Choose CSV File with longitudinal data of five coloum: ID, time, CEA, CA19-9, CA125",
                       multiple = FALSE, width = '100%',
                       accept = c(".csv"))),
      
      column(1, offset = 0.5)
    ),
    
    
    fluidRow(
      column(1, offset = 0.5),
      
      column(9, 
             tableOutput("contents")),
      
      column(1, offset = 0.5)
    ),
    
    
    fluidRow(
      column(1, offset = 0.5),
      
      column(9, align = "c",
             plotOutput("end_prediction")),
      
      column(1, offset = 0.5)
      
    ),
    
    fluidRow(
      plotOutput("dynamic_prediction"),
    )
  )
)


# Define server logic to read selected file ----
server <- function(input, output) {
  
  #导入并输出纵向数据
  dataset <- eventReactive(input$file1,{
    infile <- input$file1
    if(!is.null(infile))
      read.csv(infile$datapath, header = TRUE)
  })
  
  output$contents <- renderTable({
    return(dataset())
  },width = "100%")
  
  
  #数据处理
  f_sex <- eventReactive(input$XB,{
    ifelse(input$XB == "Male", "M", "F")
  })
  
  f_sur <- eventReactive(input$SSLJ,{
    ifelse(input$SSLJ == "Laparoscopic resection", "FQJ", "KF")
  })
  
  
  f_site <- eventReactive(input$YFBW,{
    ifelse(input$YFBW == "Colon", "JC", "ZC")
  })
  
  f_diff <- eventReactive(input$FHCD,{
    ifelse(input$FHCD == "Well", "GFH",
           ifelse(input$FHCD == "Poor-undifferentiated", "DFH",
                  ifelse(input$FHCD == "Moderate", "ZFH", "WZ")))
  })
  
  f_stage <- eventReactive(input$BLFQ,{
    ifelse(input$BLFQ == "I", "I",
           ifelse(input$BLFQ == "II", "II", "III"))
  })
  
  
  f_lny <- eventReactive(input$LYN,{
    ifelse(input$LYN == "No", "DY12", "XY12")
  })
  
  f_ac <- eventReactive(input$FZZL,{
    ifelse(input$FZZL == "Yes", "1", "0")
  })
  
  f_type <- eventReactive(input$BLLX,{
    ifelse(input$BLLX == "Yes", "NYXA", "FNYXA")
  })
  
  f_li <- eventReactive(input$MGAS,{
    ifelse(input$MGAS == "Yes", "1", "0")
  })
  
  f_pi <- eventReactive(input$SJSQF,{
    ifelse(input$SJSQF  == "Yes", "1", "0")
  })
  
  
  #数据形式转换，保证和模型内变量相同
  
  f_surv <- reactive(data.frame(ID = input$ID,  NL = as.numeric(input$NL), 
                                XB = factor(f_sex(), levels = c("M", "F")),
                                SSLJ = factor(f_sur(), levels = c("FQJ","KF")),
                                YFBW = factor(f_site(), levels = c("JC","ZC")), 
                                FHCD = factor(f_diff(), levels = c("GFH", "DFH", "WZ", "ZFH")), 
                                BLFQ = factor(f_stage(), levels = c("I", "III", "II")),
                                LYN = factor(f_lny(), levels = c("DY12","XY12")), 
                                FZZL = factor(f_ac(), levels = c("0","1")),
                                BLLX = factor(f_type(), levels = c("FNYXA", "NYXA")), 
                                MGAS = factor(f_li(), levels = c("0","1")), 
                                SJSQF = factor(f_pi(), levels = c("0","1")), 
                                JX_CEA = as.numeric(input$JX_CEA),
                                JX_CA199 = as.numeric(input$JX_CA199),
                                JX_CA125 = as.numeric(input$JX_CA125)))
  
  
  #定义预测函数
  predict <- eventReactive(input$file1,{
    infile <- input$file1
    if(!is.null(infile))
      prediction(surv = f_surv(), data_example = read.csv(infile$datapath, header = TRUE))
  })
  
  #根据定义的预测函数估计风险并画图  
  output$dynamic_prediction <- renderPlot({
    return(plot(long_example = predict()[[1]], pred_marker = predict()[[2]], pred_prob = predict()[[3]]))
  })
  
  
  output$end_prediction <- renderPlot({
    pred_plot_end <- plot_end(long_example = predict()[[1]], pred_marker = predict()[[2]], pred_prob = predict()[[3]])
    return(pred_plot_end)
  },width = 700, height = 500)
  
}

# Create Shiny app ----
shinyApp(ui, server)
