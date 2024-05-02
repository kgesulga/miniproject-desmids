plot_relabund <- function(timestep_day, comm, palette){
  fig <- dat4plots %>% 
    group_by(day_txt, community, media, temp, species, rel_abund) %>% 
    distinct(rel_abund, .keep_all = T) %>% 
    ungroup() %>% 
    filter(day %in% timestep_day,
           community %in% comm) %>%
    # plot
    ggplot(aes(fill=species_txt, y=rel_abund, x=as.factor(temp))) + 
    geom_bar(colour="black", position="stack", stat="identity") +
    facet_wrap(~media_txt, nrow = 2, ncol = 3) +
    scale_fill_manual(values = cal_palette(palette)) + #sierra1; seagrass
    ylab("Relative abundance [%]") +
    xlab("Temperature [°C]") +
    labs(fill = "Species",
         # change here the community with func argument
         title = paste("Day", timestep_day),
         subtitle = paste(comm, "community")) + 
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    # Customize facet labels
    theme(strip.background = element_rect(fill="#E4DECE"),
          strip.text.x = element_text(size = 10.5))
  
  
  # print(fig)
  
}