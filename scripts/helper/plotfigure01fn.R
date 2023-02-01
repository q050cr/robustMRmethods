

library(tibble)
library(ggplot2)

plotfigure01 <- function(data, effect_size="all", maf="all", n){
  ## data: my_data_harm
  # effect size: "all", "strong" or "weak"
  # maf: "all", "rare" or "common
  df_cum_expl_var <- tibble(n=n, 
                            cum_expl_varX = rep(NA, length(n)), 
                            cum_expl_varX2 = rep(NA, length(n)), 
                            cum_expl_varU = rep(NA, length(n)), 
                            cum_expl_varY = rep(NA, length(n)),
                            count_IVs = rep(NA, length(n))
  )
  # create maf col
  data <- data %>% 
    mutate(maf.bmi = ifelse(eaf.bmi >0.5, 1-eaf.bmi, eaf.bmi))
  
  ## filter for quadrant
  if (effect_size=="strong" & maf == "rare") {
    data <- data %>% 
      filter(abs(beta.bmi) > median(abs(beta.bmi)) &  maf.bmi < median(maf.bmi))
  } else if (effect_size=="strong" & maf == "common") {
    data <- data %>% 
      filter(abs(beta.bmi) > median(abs(beta.bmi)) &  maf.bmi > median(maf.bmi))
  } else if (effect_size=="weak" & maf == "common") {
    data <- data %>% 
      filter(abs(beta.bmi) < median(abs(beta.bmi)) &  maf.bmi > median(maf.bmi))
  } else if (effect_size=="weak" & maf == "rare") {
    data <- data %>% 
      filter(abs(beta.bmi) < median(abs(beta.bmi)) &  maf.bmi < median(maf.bmi))
  }
  
  # create data for plot
  for (i in 1:length(colnames_pval) ) {
    threshold <- 5*10^(-6)
    select_col <- colnames_pval[i]
    # calculate cum explained variance and count of IVs
    summarydat <- data %>% 
      filter(!!sym(select_col) < threshold) %>% 
      select(explained_variance_BMI, explained_variance_BMI2, explained_variance_GSD, explained_variance_GBC) %>% 
      summarise(expl_variances =across(.cols = everything(), function(x) sum(x)),
                count_IVs = n())
    # fill up df
    df_cum_expl_var[i, 2] <- summarydat$expl_variances[1]
    df_cum_expl_var[i, 3] <- summarydat$expl_variances[2]
    df_cum_expl_var[i, 4] <- summarydat$expl_variances[3]
    df_cum_expl_var[i, 5] <- summarydat$expl_variances[4]
    df_cum_expl_var[i, 6] <- summarydat$count_IVs
  }
  
  ## plot fig 01 
  scale_factor <- 4.5/284
  df_cum_expl_var <- df_cum_expl_var %>% 
    mutate(
      cum_expl_varX = cum_expl_varX *100, # in percent 
      cum_expl_varX2 = cum_expl_varX2 *100,
      cum_expl_varU = cum_expl_varU *100,
      cum_expl_varY = cum_expl_varY *100,
      count_IVs_scaled = count_IVs*scale_factor)
  
  plot_data1 <- df_cum_expl_var %>% 
    pivot_longer(cols = c(cum_expl_varX, cum_expl_varU, cum_expl_varY)) %>% 
    mutate(name = factor(name, labels = c("% var(U)", "% var(X)", "% var(Y)")),
           name = factor(name, levels = c("% var(X)", "% var(Y)", "% var(U)") )
    )
  
  N_annotation <- c(10000, 100000, 200000, 300000)
  ## plot
  plotsubtitle <- paste0("Relatively ", maf ," variants\nwith relatively ", effect_size, " effects")
  
  if( effect_size=="all" & maf=="all"  ){
    ggplot(plot_data1, aes(x=n, y=value, col=name))+
      geom_point(alpha=0.6)+
      geom_line(alpha=0.6)+
      geom_line(aes(x=n, y=count_IVs_scaled, col="Number of IVs"), linetype="dotdash")+
      scale_x_continuous(name = paste0("Size of sample ", as.roman(1)),  # 2-sample-MR or label "Study Size"
                         labels = comma,  # {scales}
                         breaks = N_annotation, #c( 10000, n[7:length(n)]),
                         limits = c(10000, 300000))+
      #ylab("% variance explained by IVs")+
      labs(title = "Figure 1", subtitle = plotsubtitle,
           col='', 
           caption = "Relationship between sample size, average number of instrumental variables (IVs) and variance of traits   explained by the IVs"
      ) +
      scale_y_continuous(
        # Features of the first axis
        name = "% variance explained by IVs",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~./scale_factor, name="Number of IVs"),
        ##limits = c(0,10)
      ) + 
      theme_few()+
      #scale_colour_few() +
      #scale_color_manual(name = "", values = c("% var(U)" = "red", 
      #                                         "% var(X)" = "blue", 
      #                                         "% var(Y)" = "darkblue", 
      #                                         "Number of IVs" = "black"
      #                                         ))
      #
      scale_color_manual(name = "", values = c("% var(U)" = ggthemes_data$few$colors$Dark[2,2][[1]], 
                                               "% var(X)" = ggthemes_data$few$colors$Dark[3,2][[1]], 
                                               "% var(Y)" = ggthemes_data$few$colors$Dark[4,2][[1]], 
                                               "Number of IVs" = "black"),
                         # order
                         breaks = c( "Number of IVs", "% var(X)", "% var(Y)", "% var(U)")
      ) -> figure01.quadrant
    
  } else{
    ggplot(plot_data1, aes(x=n, y=value, col=name))+
      geom_point(alpha=0.6)+
      geom_line(alpha=0.6)+
      geom_line(aes(x=n, y=count_IVs_scaled, col="Number of IVs"), linetype="dotdash")+
      scale_x_continuous(name = paste0("Size of sample ", as.roman(1)),  # 2-sample-MR or label "Study Size"
                         labels = comma, 
                         breaks = N_annotation, #c( 10000, n[7:length(n)]),
                         limits = c(10000, 300000))+
      #ylab("% variance explained by IVs")+
      labs(title = "Figure 1", subtitle = plotsubtitle,
           col='', 
           caption = "Relationship between sample size, average number of instrumental variables (IVs) and variance of traits   explained by the IVs"
      ) +
      scale_y_continuous(
        # Features of the first axis
        name = "% variance explained by IVs",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~./scale_factor, name="Number of IVs"),
        limits = c(0,2)
      ) + 
      theme_few()+
      #scale_colour_few() +
      #scale_color_manual(name = "", values = c("% var(U)" = "red", 
      #                                         "% var(X)" = "blue", 
      #                                         "% var(Y)" = "darkblue", 
      #                                         "Number of IVs" = "black"
      #                                         ))
      #
      scale_color_manual(name = "", values = c("% var(U)" = ggthemes_data$few$colors$Dark[2,2][[1]], 
                                               "% var(X)" = ggthemes_data$few$colors$Dark[3,2][[1]], 
                                               "% var(Y)" = ggthemes_data$few$colors$Dark[4,2][[1]], 
                                               "Number of IVs" = "black"),
                         # order
                         breaks = c( "Number of IVs", "% var(X)", "% var(Y)", "% var(U)")
      ) -> figure01.quadrant
  }
  
  print(figure01.quadrant)
  
  filenamefig01_quadrant <- paste0("../output/plots/", Sys.Date(), "_figure01_", 
                                   "EFFECT_SIZE_", effect_size,
                                   "_VARIANTS_", maf,
                                   ".png")
  ggsave(filename = filenamefig01_quadrant, plot = figure01.quadrant, 
         width = 10, height = 6, 
         units = "in"  # default
  )
  
  # store plot in global environment
  plotname <<- paste0("figure01_EFFECT_SIZE_", effect_size, "_VARIANTS_", maf)
  assign(plotname, figure01.quadrant, envir = .GlobalEnv)
  #list2env(assign(plotname, figure01.quadrant), envir = .GlobalEnv)
}

