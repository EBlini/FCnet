# #working directory - for debugging purposes
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#library(shiny)

FCnet_ui <- shiny::fluidPage(

    # Application title
    shiny::titlePanel("FCnet"),

    # Sidebar
    shiny::sidebarLayout(
        shiny::sidebarPanel(

            #x is fc matrices or volumes
            shiny::fileInput("x",
                             "Choose FC matrices or brain volumes",
                             multiple= T),

            #y is behavioral scores
            shiny::fileInput("y",
                             "Choose data to predict (one .csv file)",
                             accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            #multiple scores uploaded
            shiny::textInput("y_n",
                             "Choose data column",
                             value = 1),

            #feature reduction
            #method
            shiny::selectInput("FRmethod",
                               "Choose feature reduction method",
                               c("PCA"= "PCA",
                                 "ICA"= "ICA"),
                               selected = "PCA"),
            #number of components
            shiny::textInput("FRcomps",
                             "How many components/explained variance?",
                             value = 0.95),

            #modelling
            shiny::selectInput("alpha",
                               "Choose Regression Type",
                               c("Ridge"= 0,
                                 "Lasso"= 1,
                                 "Elastic Net"= 999),
                               selected = "Ridge"),
            shiny::fluidRow(
                shiny::column(width = 8,
                              shiny::sliderInput("FeatRange", "Features Range:",
                                   min= 10, max= 150,
                                   value= c(10, 15), step= 1)),
                shiny::column(width = 4,
                              shiny::sliderInput("FeatStep", "Features Step:",
                                   min= 1, max= 10,
                                   value= 5, step= 1))),

            shiny::selectInput("whattobp", "Back-Projection:",
                               c("Feature(s)"= 1,
                                 "Coefficient(s)"= 2),
                               selected = "Feature(s)"),

            shiny::selectInput("ctobp", "Coefficients to back-project",
                               as.character(1:150), selected = "1",
                               multiple= T),

            shiny::fluidRow(
                shiny::column(4,
                              shiny::textInput("xc", "X")),
                shiny::column(4,
                              shiny::textInput("yc", "Y")),
                shiny::column(4,
                              shiny::textInput("zc", "Z"))
            ),

            #refresh button
            shiny::submitButton(text= "Refresh"),
            shiny::downloadButton("DownloadScript", "Download Script")
        )


        ,

        # Show results
        shiny::mainPanel(

            shiny::verbatimTextOutput("Header"),
            #plotmeanFC
            shiny::plotOutput("MeanFC", width = 300, height = 300),

            #modelling + plots
            shiny::tableOutput("ModelT"),
            shiny::fluidRow(
                shiny::column(6,
                              shiny::plotOutput("ModelP1",
                                                width = 300,
                                                height = 300)),
                shiny::column(6,
                              shiny::plotOutput("ModelP2",
                                                width = 300,
                                                height = 300))),

            shiny::fluidRow(
              shiny::column(width= 12,
                            status = "primary",
                            shiny::div(style = 'overflow-x: scroll',
                                       shiny::tableOutput("Coefficients")
                            )
              )
            ),

            #backproject
            shiny::plotOutput("BPplot", width = 300, height = 300),



            shiny::verbatimTextOutput("Script")


        )

    )

)

# Define server logic required to draw a histogram
FCnet_server <- function(input, output) {

    #omit in the real app/inside package?
    #library("FCnet")


    #initial message/citations etc.
    header= paste("#Welcome!",
                  "#You are using FCnet version: ",
                   paste("#", packageVersion("FCnet")),
                  "#FCnet is available from GitHub: https://eblini.github.io/FCnet/index.html",
                  "#Bug report: elvio.blini (at) gmail.com",
                  sep= "\n")

    #render header
    output$Header= shiny::renderText({header})

    #default parameters of the app + warnings in the script
    FCnet::optionsFCnet("nested"= F)
    parallelLOO= F

    script= paste("#Generated script:\n",
                  "library('FCnet')",
                  "#Parallel processing is disabled in this app",
                  "parallelLOO= F",
                  "#Also, internal crossvalidation is not nested",
                  "optionsFCnet('nested'= F)",
                  sep= "\n")

    #1. read data

    read_data= shiny::reactive({


        shiny::req(input$x)
        shiny::req(input$y)



        #files to paste and show to the user
        yfile= input$y$name #show to user
        xfiles= input$x$name #show to user


        #internal

        y= read.csv(input$y$datapath, header= T)[,as.numeric(input$y_n)]

        read_many= function(addresses){

            if(grepl(".csv", addresses[1], fixed = T)){

                res= lapply(addresses, function(r){
                    matrix= read.csv(r, header= F)
                    matrix= as.data.frame(matrix)
                    return(matrix)

                })

            } else {



                datapathGZ= sub("gz$", "nii.gz", addresses)
                invisible(file.rename(input$x$datapath, datapathGZ))

                res= lapply(datapathGZ, function(r){



                    matrix= oro.nifti::readNIfTI(r)
                    matrix= oro.nifti::img_data(matrix)
                    return(matrix)

                })

            }

            return(res)

        }

        x= read_many(input$x$datapath)

        script_o= script
        script= paste(script_o, "\n\n",
                      "#Read files",
                      "#This bit may change according to your OS",
                      paste0("yfile= '", paste0(input$y$name, "'")),
                      paste0("y_column= ", input$y_n),
                      paste0("xfile= c('",
                             paste0(dput(input$x$name), collapse = "','"),
                             "')"),
                      "y= read.csv(yfile, header= T)[, y_column]",
                      "x= loadFC(NULL, xfile)",
                      sep= "\n"
                      )

        return(list(y= y,
                    x= x,
                    script= script))


    })

    #reduce features

    red_feat= shiny::reactive({

        shiny::req(input$x)
        shiny::req(input$y)

        #present
        RF_method= input$FRmethod

        explained_variance= as.numeric(input$FRcomps)

        RF= FCnet::reduce_featuresFC(read_data()$x,
                              RF_method,
                              Ncomp = explained_variance
                              )

        script= paste("#Reduce Features",
                      paste0("RF_method= '", RF_method, "'"),
                      paste0("explained_variance= ", explained_variance),
                      "RF= reduce_featuresFC(x, RF_method, explained_variance)",
                      sep= "\n"
        )

        return(list(RF= RF,
                    script= script))

    })

    #meanFC
    p_meanFC= shiny::reactive({

        shiny::req(input$x)
        shiny::req(input$y)

        RF= red_feat()$RF


        if(length(RF$dim)==3){

           p= FCnet::plotFC(RF$MeanFC, style = "full")

           addtoscript= "plotFC(RF$MeanFC, style = 'full')"

        } else {


            if(input$xc== "" & input$yc== "" & input$zc== ""){

                x= as.numeric(floor((RF$dim[2])/2))
                y= as.numeric(floor((RF$dim[3])/2))
                z= as.numeric(floor((RF$dim[4])/2))

            } else {
            x= as.numeric(input$xc)
            y= as.numeric(input$yc)
            z= as.numeric(input$zc)

            }

            plot.new()
            dev.control("enable")
            FCnet::plot_volume(RF$MeanFC, x= x, y= y, z= z)
            p= recordPlot()


            addtoscript= paste(
                paste0("x= ", x),
                paste0("y= ", y),
                paste0("z= ", z),
                "plot_volume(RF$MeanFC, x= x, y= y, z= z)",
                sep= "\n")
        }

        script= paste("#Plot meanFC/Volume",
                      addtoscript,
                      sep= "\n"
        )


        return(list(p= p,
                    script= script))


    })

    #plot the output meanFC
    output$MeanFC= shiny::renderPlot({p_meanFC()$p})


    #Modelling
    model_fun= shiny::reactive({

        shiny::req(input$x)
        shiny::req(input$y)

        RF= red_feat()$RF

        alpha= as.numeric(input$alpha)
        if(alpha== 999)(alpha=seq(0, 1, 0.1))

        min_comp= input$FeatRange[1]
        max_comp= input$FeatRange[2]
        if(max_comp> ncol(RF$Weights)){

            max_comp= ncol(RF$Weights)
        }
        step= input$FeatStep

        model= FCnet::FCnetLOO(read_data()$y,
                        RF,
                        alpha = alpha,
                        cv_Ncomp = seq(min_comp,
                                       max_comp,
                                       step))

        p1= FCnet::plotFCnet(model, plot_labels = F)
        p2= FCnet::plotFCnet(model, "coefficients")



        modeltres= data.frame(R2= model$R2,
                              alpha= model$alpha,
                              lambda= model$lambda,
                              Features= ncol(RF$Weights),
                              k= length(model$N_comp),
                              NonZero= sum(model$coeffs$Coefficient!=0),
                              N_obs= length(model$y))

        script= paste("#Modelling",
                      paste0("alpha= ", alpha),
                      paste0("min_comp= ", min_comp),
                      paste0("max_comp= ", max_comp),
                      paste0("step= ", step),
                      "model= FCnetLOO(y, RF, alpha = alpha, cv_Ncomp = seq(min_comp, max_comp, step))",
                      sep= "\n")

        return(list(model= model,
                    modeltres= modeltres,
                    p1= p1, p2= p2,
                    script= script))
    })

    #tableresults
    output$ModelT= shiny::renderTable({model_fun()$modeltres})
    #plot results
    output$ModelP1= shiny::renderPlot({model_fun()$p1})
    output$ModelP2= shiny::renderPlot({model_fun()$p2})

    #coefficients
    output$Coefficients=shiny::renderTable({
      reshape2::dcast(model_fun()$model$coeffs, . ~ Feature,
                      value.var="Coefficient")
        })

    #backproject
    bp_fun= shiny::reactive({

        shiny::req(input$x)
        shiny::req(input$y)

        RF= red_feat()$RF
        model= model_fun()$model

        if(input$whattobp== "1"){
        coeffs= as.numeric(input$ctobp)
        addtoscript1= paste("coeffs= c('",
                           paste(ifelse(length(dput(input$ctobp))>1,
                                        dput(input$ctobp), input$ctobp),
                           collapse = ifelse(length(dput(input$ctobp))>1,"','", "")
                           ), "')")

        } else {
            coeffs= as.numeric(input$ctobp) #yeah redundant
            #grab coefficients
            vc= rep(0, nrow(model_fun()$model$coeffs)-1) #minus intercept
            vc[coeffs]= model_fun()$model$coeffs$Coefficient[coeffs + 1]
            coeffs= vc
            addtoscript1= paste("coeffs= c('",
                               paste(dput(coeffs),
                                     collapse = "','"), "')")

        }


        BP= FCnet::backprojectFCnet(coeffs,
                             RF)



        addtoscript2= "BP= backprojectFCnet(coeffs, RF)"



        if(length(RF$dim)==3){

            p= FCnet::plotFC(BP, style = "full")

            addtoscript3= "plotFC(BP, style = 'full')"

        } else {

            if(input$xc== "" & input$yc== "" & input$zc== ""){

                x= as.numeric(floor((RF$dim[2])/2))
                y= as.numeric(floor((RF$dim[3])/2))
                z= as.numeric(floor((RF$dim[4])/2))

            } else {
                x= as.numeric(input$xc)
                y= as.numeric(input$yc)
                z= as.numeric(input$zc)

            }

            plot.new()
            dev.control("enable")
            FCnet::plot_volume(BP, x= x, y= y, z= z)
            p= recordPlot()


            addtoscript3= paste(
                paste0("x= ", x),
                paste0("y= ", y),
                paste0("z= ", z),
                "plot_volume(BP, x= x, y= y, z= z)",
                sep= "\n")
        }

        script= paste("#Back-Projection",
                      addtoscript1,
                      addtoscript2,
                      addtoscript3,
                      sep= "\n"
        )


        return(list(p= p,
                    script= script))




    })

    output$BPplot= shiny::renderPlot({bp_fun()$p})


    merge_script= reactive({

        req(input$x)
        req(input$y)

        final_script= paste(header,
                            read_data()$script,
                            red_feat()$script,
                            p_meanFC()$script,
                            model_fun()$script,
                            bp_fun()$script,
                            sep= "\n\n"
                            )

        return(final_script)

    })

    output$Script= shiny::renderText({merge_script()})

    output$DownloadScript= shiny::downloadHandler(

        filename = "FCnet_script",

        content = function(filename) {

            obj= merge_script()

            write.table(obj,
                      filename)
        }
    )





}

# Run the application
#' FCnet (experimental) shiny Graphic User Interface
#'
#' This function launches an experimental GUI for FCnet
#' implementing a slightly simplified pipeline.

#' @export

run_FCnet= function(){

    shiny::runApp(
        shiny::shinyApp(ui = FCnet_ui, server = FCnet_server)
    )

}


