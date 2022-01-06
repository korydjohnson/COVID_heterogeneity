######################################################3
# Plotting Functions for Multivariant Simulation
######################################################3

# colour_palette <-
#   c("#56B4E9",
#     "#009E73",
#     "#F0E442",
#     "#0072B2",
#     "#D55E00",
#     "#CC79A7",
#     "#E69F00")


### INCIDENCE PLOTS
plot_incidence = function(data,
                          reported_cases,
                          nu = 13,
                          variants_to_plot = c("W", "A", "B", "G", "D"),
                          plot_cutoff = 200,
                          scale_by = 10 ^ 5 / 8932664,
                          facet_by = NULL,
                          facet_rows = 1,
                          detRatio = 1/1.4) {
  data[c('I_2.5%', 'I_25%', 'I_50%', 'I_75%', 'I_97.5%')] =
    data[c('I_2.5%', 'I_25%', 'I_50%', 'I_75%', 'I_97.5%')] * scale_by * detRatio
  
  incidence_plot = data %>%
    filter(variant %in% c(variants_to_plot, "Total")) %>%
    mutate(variant = factor(variant, levels =  c(variants_to_plot, "Total"))) %>%
    ggplot(aes(x = date)) +
    #columns for reported cases
    geom_area(
      data = reported_cases,
      aes(x = day, y = cases_ra * scale_by),
      fill = "#8C8C8C",
      alpha = 0.25
    ) +
    ##RIBBONS
    geom_ribbon(aes(
      ymin = `I_25%`,
      ymax = `I_75%`,
      fill = variant
    ), alpha = .5) +
    geom_ribbon(aes(
      ymin = `I_2.5%`,
      ymax = `I_97.5%`,
      fill = variant
    ), alpha = .15) +
    ##Line for the forecast median
    geom_line(aes(
      y = `I_50%`,
      color = variant,
      size = variant,
      # alpha = variant,
      linetype = variant
    )) +
    ## define line widths
    scale_size_manual(values = c(rep(0.8, length(variants_to_plot)), 1.2)) +
    # scale_alpha_manual(values = c(rep(1, length(variants_to_plot)), 0.8)) +
    # scale_linetype_manual(values = c(rep("solid", length(variants_to_plot)), "dotted")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(),
      aspect.ratio = 0.5
    ) +
    ylab("daily incidence")  +
    coord_cartesian(ylim = c(0, plot_cutoff)) +
    theme_bw()
  
  if (!is.null(facet_by)) {
    incidence_plot = incidence_plot + facet_wrap(vars(!!sym(facet_by)), nrow = facet_rows, scales = "free_x")
  }
  
  incidence_plot
}

plot_stat = function(stat, data_stats=NULL, file_list=NULL, facet_by = NULL,
                     variants_to_plot = "", end_date = as.Date("2022-03-01")) {
  if (is.null(data_stats) & is.null(file_list)) {
    stop("Need either data object or data_list.")
  }
  if (is.null(data_stats)) {
    forecasting_date = ifelse(str_detect(file_list$variant, "Delta"), 
                              "2021-09-01", "2021-08-08")
    data_stats = do.call(read_samples, file_list) %>% 
      process_simData %>% 
      .$stats %>% 
      mutate(date = day + as.Date(forecasting_date)) %>%
      filter(date <= end_date)
  }
  if (stat == "rhoV") {
    data_stats = data_stats %>% 
      select(date, starts_with("rhoV")) %>% 
      pivot_longer(starts_with("rhoV"), names_prefix = "rhoV.") %>% 
      separate(name, into = c("variant", "q"), sep="_") %>% 
      pivot_wider(names_from="q", names_prefix = "rhoV_") %>% 
      filter(variant %in% variants_to_plot)
    plt = ggplot(data_stats, aes(x = date, fill=variant)) +
      geom_line(aes(y = !!sym(paste0(stat, "_50%")), color=variant, linetype=variant))
  } else {
    plt = ggplot(data_stats, aes(x=date)) +
      geom_line(aes(y = !!sym(paste0(stat, "_50%"))))
  }
  plt = plt +
    ##RIBBONS
    geom_ribbon(
      aes(ymin = !!sym(paste0(stat, "_25%")), ymax = !!sym(paste0(stat, "_75%"))),
      alpha = .4
    ) +
    geom_ribbon(
      aes(ymin = !!sym(paste0(stat, "_2.5%")),ymax = !!sym(paste0(stat, "_97.5%"))),
      alpha = .25
    ) +
    ##Line for the forecast median
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 0.3) +
    ylab(stat) +
    theme_bw()
  
  if (!is.null(facet_by)) {
    plt +facet_wrap(vars(!!sym(facet_by)), scales = "free_x")
  } else {
    plt
  }
}

add_vert_lines = function(plt, forecasting_date) {
  # vertical lines for starting days
  plt +
    geom_vline(
      aes(xintercept = forecasting_date-12),
      linetype = "dashed",
      alpha = 0.25
    ) +
    geom_vline(
      aes(xintercept = forecasting_date),
      linetype = "dashed"
    )
}

make_side_by_side_plot = function(data_list,
                                  cases,
                                  data_folder = "",
                                  end_date = as.Date("2022-03-01"),
                                  nu = 13,
                                  # save_plot = T,
                                  scale_by = 10 ^ 5 / 8932664) {
  variant_names = c("W", "A", "B", "G", "D", "O")
  variant_names_full = c("WT", "Alpha", "Beta", "Gamma", "Delta", "Omega")
  colour_palette = c("#E69F00",
                     "#56B4E9",
                     "#009E73",
                     "#F0E442",
                     "#0072B2",
                     "#D55E00",
                     "#CC79A7")
  # colors = ggthemes::colorblind_pal()(8)[c(3,7,4,2,1)]
  colours = set_names(colour_palette, c(variant_names_full, "Total"))
  shapes = c("twodash",
             "longdash",
             "dotdash",
             "F1",
             "solid",
             "dashed",
             "dotted")
  linetypes = set_names(shapes, c(variant_names_full, "Total"))
  
  info_right = data_list$info_right
  info_left = data_list$info_left
  
  upper_left = upper_right = data_list$upper
  lower_left = lower_right = data_list$lower
  
  #data left
  data_left = do.call(read_samples, data_list$file_left) %>% 
    process_simData
  #data right
  data_right = do.call(read_samples, data_list$file_right) %>% 
    process_simData
  
  data_incidence =
    rbind(
      data_left$incidence %>%
        mutate(info = info_left,
               date = day + data_list$forecasting_date),
      data_right$incidence %>%
        mutate(info = info_right,
               date = day + data_list$forecasting_date)
    ) %>%
    mutate(info = factor(info, levels = c(info_left, info_right))) %>%
    filter(date <= end_date)
  
  data_stats =
    rbind(
      data_left$stats %>%
        mutate(info = info_left,
               date = day + data_list$forecasting_date),
      data_right$stats %>%
        mutate(info = info_right,
               date = day + data_list$forecasting_date)
    ) %>%
    mutate(info = factor(info, levels = c(info_left, info_right))) %>%
    filter(date <= end_date)
  
  
  
  for (i in 1:length(variant_names)) {
    data_incidence[data_incidence$variant == variant_names[i], 2] =
      variant_names_full[i]
  }
  
  reported_cases_left = cases %>%
    filter(date >= data_list$forecasting_date - nu + 1) %>%
    mutate(day = as.integer(date - data_list$forecasting_date)) %>%
    select(day = date, cases_total, cases_ra)
  
  reported_cases_right = cases %>%
    filter(date >= data_list$forecasting_date - nu + 1) %>%
    mutate(day = as.integer(date - data_list$forecasting_date)) %>%
    select(day = date, cases_total, cases_ra)
  
  reported_cases =
    rbind(
      reported_cases_left %>% mutate(info = info_left),
      reported_cases_right %>% mutate(info = info_right)
    ) %>%
    mutate(info = factor(info, levels = c(info_left, info_right)))
  
  colour_palette = colour_palette_fill =
    colours[c(data_list$variants_to_plot, 'Total')]
  colour_palette_fill["Total"] = "transparent"
  
  incidence_plot  = plot_incidence(
    data = data_incidence,
    reported_cases = reported_cases,
    nu = 13,
    variants_to_plot = data_list$variants_to_plot,
    plot_cutoff = data_list$plot_cutoff,
    #cause we plot cases/10^5
    scale_by = scale_by,
    facet_by = "info",
    facet_rows = 1
  ) +
    scale_fill_manual(values = colour_palette_fill) +
    scale_colour_manual(values = colour_palette) +
    scale_linetype_manual(values = linetypes[c(data_list$variants_to_plot, 'Total')]) +
    guides(size = "none") +
    theme(legend.box = "horizontal") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    labs(fill = "Variant",
         color = "Variant",
         linetype = "Variant") +
    #background_grid() +
    theme(legend.position = "bottom",
          aspect.ratio = 0.5) +#,
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=.8)) +
    ylab(expression(paste("Daily Incidence per ", 10 ^ 5))) +
    xlab("Date") +
    # horizontal lines for strict case number restrictions
    geom_hline(
      data = filter(data_incidence, info == info_left),
      aes(yintercept = upper_left * scale_by),
      alpha = 0.5,
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_hline(
      data = filter(data_incidence, info == info_left),
      aes(yintercept = lower_left * scale_by),
      alpha = 0.5,
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_hline(
      data = filter(data_incidence, info == info_right),
      aes(yintercept = upper_right * scale_by),
      alpha = 0.5,
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_hline(
      data = filter(data_incidence, info == info_right),
      aes(yintercept = lower_right * scale_by),
      alpha = 0.5,
      size = 0.5,
      linetype = "dotted"
    )
  incidence_plot = add_vert_lines(incidence_plot, data_list$forecasting_date)
  
  
  mit_plot = plot_stat("mit", data_stats, facet_by = "info") + 
    ylim(c(0, 1)) +
    scale_x_date(date_breaks = "months", date_labels = "%b",
                 limits = range(incidence_plot$data$date)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      aspect.ratio = 0.3
    ) +
    ylab(expression(paste(tilde(M)[t], "= 1-", M[t])))
  mit_plot = add_vert_lines(mit_plot, data_list$forecasting_date)
  
  
  R_plot = plot_stat("R_e_t", data_stats, facet_by = "info") +
    scale_y_continuous(breaks = seq(0, 2, 0.2)) +
    scale_x_date(date_breaks = "months", date_labels = "%b",
                 limits = range(incidence_plot$data$date)) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      aspect.ratio = 0.3
    ) +
    ylab(expression(paste(hat(R)[e][","][t]))) +
    geom_hline(aes(yintercept = 1), alpha = 0.25) +
    #geom_hline(aes(yintercept = 0.8), alpha = 0.5,linetype = "dotted") +
    geom_hline(
      data = filter(data_stats, info == info_left),
      aes(yintercept = data_list$upper_R),
      alpha = 0.5,
      size = 0.5,
      linetype = "dotted"
    )
  R_plot = add_vert_lines(R_plot, data_list$forecasting_date)
  
  g = arrangeGrob(rbind(
    ggplotGrob(R_plot),
    ggplotGrob(mit_plot),
    ggplotGrob(incidence_plot)
  ))
  ggsave(data_list$plot_file, plot = g, width = 7, height = 5.6,  dpi=600)
}

make_prevalence_plots = function(data_list, data_folder = "",
                                 variant_data = read.csv("data/ages_variant_data2.csv")) {
  KW_dates = as.Date("2021-01-10") + (0:51) * 7
  variant_data %<>% 
    mutate(date = KW_dates[KW])
  
  for (file_list in list(data_list$file_left, data_list$file_right)) {
    plotFile = str_replace(data_list$plot_file, "\\.", 
                           paste0("_", file_list$control, "\\."))
    
    prevalence_data = do.call(read_samples, file_list) %>% 
      process_simData %>%
      .$incidence %>% 
      mutate(date = day + data_list$forecasting_date) %>%
      filter(variant == "D")
    
    prevalence_plot = prevalence_data %>%
      filter(date <= as.Date("2021-08-08"), date >= data_list$forecasting_date - 7) %>%
      ggplot(aes(x = date)) +
      ##RIBBONS
      geom_ribbon(
        aes(ymin = `prevalence_25%`, ymax = `prevalence_75%`),
        fill = "#0072B2",
        alpha = .25
      ) +
      geom_ribbon(
        aes(ymin = `prevalence_2.5%`, ymax = `prevalence_97.5%`),
        fill = "#0072B2",
        alpha = .1
      ) +
      ##Line for the forecast median
      geom_line(aes(y = `prevalence_50%`, colour = "Simulation"),
                alpha = 0.8,
                size = 1) +
      geom_point(
        data = variant_data %>%
          filter(KW >= 23, KW <= 31),
        aes(x = date, y = D, size = "Observed Proportion (Weekly)"),
        colour = "#CC79A7"
      ) +
      geom_vline(aes(xintercept = data_list$forecasting_date), linetype = "dashed") +
      theme_bw() +
      theme(aspect.ratio = 0.3, legend.position = "bottom") +
      ylab("Delta Prevalence") +
      xlab("Date") +
      guides(fill = "none") +
      labs(size = "", colour = "") +
      scale_colour_manual(values = "#0072B2") +
      scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %y") +
      scale_y_continuous(breaks = seq(0, 1, by = 0.2), lim = c(0, 1))
    
    ggsave(plotFile, plot = prevalence_plot, width = 7, height = 3,  dpi=600)
  }
}

make_rho_plot = function(data_list, toReturn=F) {
  colour_palette = set_names(
    c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
    c("W", "A", "B", "G", "D", "O", "T")
  )
  shapes = set_names(
    c("twodash", "longdash", "dotdash", "F1", "solid", "dashed", "dotted"),
    c("W", "A", "B", "G", "D", "O", "T")
  )
  
  variant_names = c("A", "B", "G", "D", "O")
  variant_names_full = c("Alpha", "Beta", "Gamma", "Delta", "Omega")
  plotList = list()
  
  for (file_list in list(data_list$file_left, data_list$file_right)) {
    plotFile = str_replace(data_list$plot_file, "\\.", 
                           paste0("_", file_list$control, "\\."))
    
    plt = plot_stat("rhoV", file_list = file_list, 
                    variants_to_plot = variant_names) +
      theme(legend.position = "bottom",
            aspect.ratio = 0.5) +
      scale_colour_manual(values = colour_palette[variant_names],
                          labels = variant_names_full) +
      scale_fill_manual(values = colour_palette[variant_names],
                        labels = variant_names_full) +
      scale_linetype_manual(values = shapes[variant_names],
                            labels = variant_names_full) +
      xlab("Date") +
      scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
      ylab(expression(paste(rho[t] ^ V))) +
      labs(colour = "Variant", fill = "Variant", linetype="Variant")
    if (!toReturn) {
      ggsave(plotFile, plot = plt, width = 5, height = 5*2/3,  dpi=600)
    }
    plotList[[length(plotList)+1]] = plt
  }
  if (toReturn) {
    plotList
  }
}

print_difference = function(compareTib, colName, data_list) {
  reactive = compareTib %>% 
    filter(control=="reactive") %>% 
    pull(!!sym(colName))
  proactive = compareTib %>% 
    filter(control=="proactive") %>% 
    pull(!!sym(colName))
  testOut = wilcox.test(reactive, proactive, conf.int = TRUE)
  
  print(unlist(data_list$file_left))
  print(colName)
  print(testOut)
  
  print("Proactive PI")
  print(quantile(proactive, probs = c(.025, .5, .975)))
  print("Reactive PI")
  print(quantile(reactive, probs = c(.025, .5, .975)))
  
  # jmuOutlier::quantileCI(x=x, probs=p, conf.level=0.95)[1,c("lower","upper")])
  # jmuOutlier::quantileCI
}

make_comparison_plots = function(data_list, toReturn=F) {
  compareTib = list(do.call(read_samples, data_list$file_left),
                    do.call(read_samples, data_list$file_right)) %>% 
    unlist(recursive=F) %>% 
    compare_simData
  simLength = length(compareTib$mit[[1]]) - 1
  
  infPlot = ggplot(compareTib) +
    geom_boxplot(aes(control, infections)) +
    # facet_grid(rows = vars(bounds), scales="free_y") +
    ylab(expression("Total Infections per 10"^5)) +
    xlab("Controller Type") +
    # ggtitle("Delta Simulation Starting 2021-06-12") +
    theme_bw()
  
  print_difference(compareTib, "infections", data_list)
  if (!toReturn) {
    plotFile = str_replace(data_list$plot_file, "\\.", "_Infections\\.")
    ggsave(plotFile, plot = infPlot, width = 4, height = 4,  dpi=600)
  }
  
  quarantinePlot = ggplot(compareTib) +
    geom_boxplot(aes(control, quarantine)) +
    # facet_grid(rows = vars(bounds), scales="free_y") +
    ylab(expression("Peak Quarantined Cases per 10"^5)) +
    xlab("Controller Type") +
    # ggtitle("Delta Simulation Starting 2021-06-12") +
    theme_bw()
  
  print_difference(compareTib, "quarantine", data_list)
  if (!toReturn) {
    plotFile = str_replace(data_list$plot_file, "\\.", "_Quarantine\\.")
    ggsave(plotFile, plot = quarantinePlot, width = 4, height = 4,  dpi=600)
  }
  
  aveMit = compareTib %>% 
    group_by(control) %>% 
    summarise(mit = list(matrix(unlist(mtilde), nrow=n(), byrow=T) %>% colMeans()),
              mstar = list(matrix(unlist(mstar), nrow=n(), byrow=T) %>% colMeans())) %>% 
    unnest(everything())
  
  mitPlot = aveMit %>% 
    mutate(difference = mstar-mit) %>% 
    ggplot(aes(difference, color=control, fill=control)) +
    geom_density(alpha=.3) +
    geom_vline(aes(xintercept=0)) +
    scale_color_colorblind() +
    scale_fill_colorblind() +
    theme_bw() +
    theme(legend.position = c(.85,.85)) +
    xlab(expression(paste(tilde(M)[t]^"*", " - ", tilde(M)[t]))) +
    labs(fill = "Controller", color = "Controller")
  
  
  if (!toReturn) {
    plotFile = str_replace(data_list$plot_file, "\\.", "_Mitigation\\.")
    ggsave(plotFile, plot = mitPlot, width = 4, height = 4,  dpi=600)
  }
  if(toReturn) {
    list(infPlot=infPlot, quarantinePlot=quarantinePlot, mitPlot=mitPlot)
  }
}

mergePlots = function(..., captions=NULL, fileName = "test", height=4, keepScale=F) {
  plotList = list(...)
  library(cowplot)
  if (is.null(captions)) {
    captions = rep("", length(plotList))
  }
  g = plotList[[2]] +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.justification="center")
  legend = get_legend(g)
  if (keepScale) {
    rng = list()
    for (i in seq_along(plotList)) {
      rng[[i]] = layer_scales(plotList[[i]])$y$range$range
    }
    rngFull = range(unlist(rng))
    for (i in seq_along(plotList)) {
      plotList[[i]] = plotList[[i]] +
        ylim(rngFull)
    }
  }
  for (i in seq_along(plotList)) {
    plotList[[i]] = plotList[[i]] +
      theme(legend.position = "none",
            plot.caption = element_text(hjust = 0)) +
      labs(caption = paste0("(", letters[i], ") ", captions[i]))
  }
  plotFile = paste0("plots/", fileName)
  grid1 = plot_grid(plotlist=plotList, nrow = 1, axis="tb", align="v")
  plt = plot_grid(grid1, legend, ncol=1, rel_heights = c(1, .1))#, scale=c(1,.5))
  ggsave(plotFile, plot = plt, width = 4*length(plotList), 
         height = height,  dpi=600)
}
