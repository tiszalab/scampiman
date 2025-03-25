#!/usr/bin/env Rscript

### script to make a plot of amplicon coverage for samples from 1 or more projects.
### It takes 2 arguments: 
### 1) directory to find scampiman output files and index sheets
### Projects are deciphered by subdirectories within arg[1] if any.
### Index sheets have to be placed into directories with its samples
### Index columns: Sample ID,Barcode ID,Primer,Forward Sequence,Reverse Sequence
### 2) name of file to save image of plots. Should have .png or .pdf extension



args = commandArgs(trailingOnly=TRUE)

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(rprojroot))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(readxl))


# load tables
## find files

table_files <- list.files(
  args[1], 
  recursive = T,
  pattern = "*amplicontable.tsv", 
  full.names = TRUE)

## make big dt
if (exists("amp_dataset")){
  rm(amp_dataset)
}
for (tabf in table_files){
  basef <- basename(tabf)
  projID <- basename(dirname(tabf))
  
  if (!exists("amp_dataset")){
    amp_dataset <- fread(
      tabf,
      sep = "\t",
      header = T
    ) %>%
      mutate(project_ID = projID)
  }
  else if (exists("amp_dataset")){
    tmp_dataset <- fread(
      tabf,
      sep = "\t",
      header = T
    ) %>%
      mutate(project_ID = projID)
    
    amp_dataset <- rbind(
      amp_dataset,
      tmp_dataset
    )
    rm(tmp_dataset)
  }
  
}

# load index sheets

index_files <- list.files(
  args[1], 
  recursive = T,
  pattern = "*.xlsx", 
  full.names = TRUE)

#index_files

## make big dt
if (exists("index_dt")){
  rm(index_dt)
}
for (indexf in index_files){
  basef <- basename(indexf)
  projID <- basename(dirname(indexf))
  
  if (!exists("index_dt")){
    index_dt <- read_excel(
      indexf
    ) %>%
      mutate(project_ID = projID,
             sample_ID = gsub("NB", "barcode",`Barcode ID`))
  }
  else if (exists("index_dt")){
    tmp_dataset <- read_excel(
      indexf
    ) %>%
      mutate(project_ID = projID,
             sample_ID = gsub("NB", "barcode",`Barcode ID`))
    
    index_dt <- rbindlist(
      list(index_dt,
      tmp_dataset),
      fill = TRUE
    ) 
    rm(tmp_dataset)
  }
  
}

#index_dt

amp_index_dt <- merge(
  amp_dataset,
  index_dt,
  by = c("sample_ID", "project_ID")
)

heatp <- amp_index_dt %>%
  ggplot(aes(
    x = amplicon_number,
    y = `Sample ID`,
    fill = log10(amplicon_reads)
  )
  ) +
  geom_tile(color = NA) +
  facet_grid(vars(project_ID),
             vars(accession),
             scales = "free",
             space = "free"
  ) +
  scale_fill_gradient2(
    low = "gold",
    mid = "#54D8B1",
    high = "#175149",
    na.value = "grey90",
    midpoint = 2
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


barp <- amp_index_dt %>%
  group_by(
    project_ID, `Sample ID`, accession
  ) %>%
  summarize(
    amps_covered = sum(amplicon_reads  != 0),
    avg_cov = mean(amplicon_reads),
    n_amps = n()
  ) %>%
  mutate(perc_amps = amps_covered/n_amps) %>%
  ggplot(aes(
    x = perc_amps,
    y = `Sample ID`,
    fill = avg_cov
  )
  ) +
  geom_col(color = NA) +
  geom_text(aes(label = scales::percent(perc_amps,
                                        accuracy = 1),
                x = 0),
            nudge_x = 0.1,
            size = 3,
            hjust = 0) +
  scale_fill_gradient(
    low = "#90D4CC",
    high = "#BD3027"
  ) +
  facet_grid(vars(project_ID),
             vars(accession),
             scales = "free"
  ) +
  theme_bw() +
  theme(axis.text.y.left = element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "bottom")


combp <- cowplot::plot_grid(heatp, 
                            barp, 
                            align = "h", 
                            nrow = 1, 
                            rel_widths = c(70/100, 30/100)
                            )

dir.create(file.path(args[1], "charts"), showWarnings = FALSE)

ggsave(
  combp,
  file = sprintf(
    "%s/charts/%s",
    args[1],
    args[2]
  ),
  width = 15,
  height = 12
)






