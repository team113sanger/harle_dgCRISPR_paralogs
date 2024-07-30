### This script was used for image analysis 
### Author: Victoria Harle

# Purpose:
# Determine and remove outliers
# Robust Z-transform data with respect to controls

# Load dependencies
library(tidyverse)
library(ggrepel)
library(scales) 
library(ggpubr)
library(egg)
library(pheatmap)
library(plotly)
library(gridExtra)

# Top level directory - assumes you are in the directory with the imaging data and scripts
top_dir <- getwd()

# Source helper functions
source(file.path(top_dir, 'SCRIPTS', 'helper.R'))

# Read in raw data list from RDS
data.list <- readRDS(file = file.path(top_dir, 'RESULTS', 'raw_processed_plate_list_with_medians.rds'))

# Collate raw data set into a single data frame ---------------------------
# Collate nested list into single data frame
# Remove blank wells
data.df <- convert_nested_list_to_df(data.list)

# Remove wells with analysed fields less than 10 --------------------------
# Remove wells which have less than 10 analysed fields
data.df <- data.df |> filter(`Number of Analyzed Fields` >= 10)

# Remove wells which have less than 1000 objects --------------------------
# Remove wells which have less than 1000 objects
data.df <- data.df |> filter(`Cells Final - Number of Objects` >= 1000)

#############REMOVE PROLIF UNSTAINTED OUTLIERS - several low percentages stick out from two plates only as poorly stained (outliers from group so removed)
data.df <- data.df |> filter(`% Proliferative` > 5)

############Exclude the PCA outlier control wells for later scaling
### Question - do we need to include the script to show these wells should be excluded? - Victoria??
data.df <- data.df |> filter(!c(Replicate =="N1" & Plate =="Plate_4" & (Well == '2,10'| Well == '2,11'))) |> 
  filter(!c(Replicate =="N2" & Plate =="Plate_4" & (Well == '7,10'| Well == '7,11')))

# Scale wells by replicate contribution per plate -------------------------
# Get total number of objects per well for all replicates - Note this is scaling one, for the objects (no. of cells). Think normalisation by cell number so the proportions are comparable across the whole plate.
# Calculate scaling factor as proportion contribution of each replicate to total
object_summary <- data.df |> 
  select(Plate:Replicate, `Cells Final - Number of Objects`) |>
  group_by(Plate, Replicate) |>
  summarise('total_objects_per_well' = sum(`Cells Final - Number of Objects`), .groups = 'keep') |>
  group_by(Plate) |>
  mutate('total_objects' = sum(.data$total_objects_per_well)) |>
  mutate('scaling_factor' = total_objects_per_well / total_objects)

# Divide each well value by the corresponding scaling factor (between-plate scaling)
# Note: there is no within-plate scaling to controls here
scaled_objects.narrow <- data.df |>
  select(Plate:Replicate, `Non-Proliferative - Number of Objects`:`Cells Final - Number of Objects`) |>
  pivot_longer(cols = `Non-Proliferative - Number of Objects`:`Cells Final - Number of Objects`, values_to = 'raw') |>
  left_join(object_summary |> select(Plate, Replicate, scaling_factor), by = c('Plate', 'Replicate'), multiple = "all") |>
  mutate('scaled' = raw / scaling_factor)

# Spread data set - This is the data to use for plotting percentages of cells - however remember to remove the unclassified and plot only the classified populations
scaled_objects.wide <- scaled_objects.narrow |>
  pivot_wider(names_from = name, values_from = c(raw, scaled)) 

# Remove genes 
scaled_objects.wide <- scaled_objects.wide |> 
  filter(!Group_Target %in% c('TTC7A_TTC7B', 'INTS6_INTS6L', 'EAF1_EAF2'))

#Calculate sum of classified objects only (i.e. remove those unclassified to allow us to work out percentage of classified cells)
## Question - should this be done before scaling, i.e. does removing the unclassified objects affect the scaling at all?? - Victoria?? 
Classified <- scaled_objects.wide %>% rowwise() %>% 
  mutate(total_scaled_classified_objects = (sum(`scaled_Non-Proliferative - Number of Objects`,
                                                `scaled_Proliferative - Number of Objects`,
                                                `scaled_Apoptotic - Number of Objects`,
                                                `scaled_Enlarged - Number of Objects`)))

#Work out classified percentage populations 
Average_percentage_scaled <- Classified  %>% rowwise() %>% 
  mutate(Percentage_Non_proliferative = (sum((`scaled_Non-Proliferative - Number of Objects`)/ sum(total_scaled_classified_objects))*100),                                                        
         Percentage_Proliferative = (sum((`scaled_Proliferative - Number of Objects`)/ sum(total_scaled_classified_objects))*100),                                               
         Percentage_Apoptotic = (sum((`scaled_Apoptotic - Number of Objects`)/ sum(total_scaled_classified_objects))*100),                                               
         Percentage_Enlarged = (sum((`scaled_Enlarged - Number of Objects`)/ sum(total_scaled_classified_objects))*100)) 

#List targets in order to plot
names<-c("Control_1", "Control_2", "Control_1|Control_2",  "CNOT7","CNOT8","CNOT7|CNOT8","CCNL1","CCNL2","CCNL1|CCNL2",
         "ASF1A","ASF1B","ASF1A|ASF1B", "SLC25A37","SLC25A28","SLC25A37|SLC25A28", "GDI1","GDI2","GDI1|GDI2", 
         "PDS5A","PDS5B","PDS5A|PDS5B","SAR1A","SAR1B","SAR1A|SAR1B" ,"SEC23A","SEC23B","SEC23A|SEC23B",
         "EAF1","EAF2","EAF1|EAF2" ,"INTS6","INTS6L","INTS6|INTS6L","TTC7A","TTC7B","TTC7A|TTC7B", 'Parental')

Plot_order<-Average_percentage_scaled %>% arrange(factor(Target, levels=names))


#### Separate plots of percentage populations classified, scaled
Plot_percentage_non_prolif <- ggplot(Plot_order, aes(x=fct_inorder(Target), y=Percentage_Non_proliferative, fill=Group_Target)) + 
  geom_boxplot() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') + 
  labs(title='Percentage Not Proliferative', y='Percentage Not Proliferative', x='Target')+ylim(0,100)

Plot_percentage_prolif <- ggplot(Plot_order, aes(x=fct_inorder(Target), y=Percentage_Proliferative, fill=Group_Target)) + 
  geom_boxplot() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') + 
  labs(title='Percentage Proliferative Cells', y='Percentage Proliferative', x='Target')+ylim(0,100)

Plot_percentage_apoptotic <- ggplot(Plot_order, aes(x=fct_inorder(Target), y=Percentage_Apoptotic, fill=Group_Target)) + 
  geom_boxplot() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') + 
  labs(title='Percentage Apoptotic Cells', y='Percentage Apoptotic', x='Target') +ylim(0,100)

Plot_percentage_enlarged <- ggplot(Plot_order, aes(x=fct_inorder(Target), y=Percentage_Enlarged, fill=Group_Target)) + 
  geom_boxplot() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') + 
  labs(title='Percentage Enlarged Cells', y='Percentage Enlarged', x='Target') +ylim(0,15)

#Arrange percentage proliferation plots 
combined_percentage_classified_figure <- grid.arrange(Plot_percentage_non_prolif,
                                                      Plot_percentage_prolif,
                                                      Plot_percentage_apoptotic,
                                                      Plot_percentage_enlarged, nrow=2)

# Save combined plot to file ----------------------------------------------
ggsave(file.path(top_dir, 'PLOTS', 'combined_percentage_classified_figure.png'), combined_percentage_classified_figure , dpi = 300, width = 12, height = 10)

#######Percentage hits stacked bar plot generation for figure -------------------------------------------------------------------------------------
######Get scaled data and average per target
### Remove Parental cells
Average_percentage_scaled <- Average_percentage_scaled %>% filter(!Target == 'Parental')

colnames(Average_percentage_scaled)
Percentage_scaled<- Average_percentage_scaled %>% group_by(Target) %>% mutate(`Non Proliferative` = mean(Percentage_Non_proliferative),
                                                                              `Proliferative` = mean(Percentage_Proliferative),
                                                                              `Apoptotic` = mean(Percentage_Apoptotic),
                                                                              `Enlarged` = mean(Percentage_Enlarged)) %>%
  select(c("Target", "Group_Target", "Non Proliferative" , "Proliferative", "Apoptotic", "Enlarged")) %>% unique()

Percentage_scaled_pivot <- Percentage_scaled %>% pivot_longer(c(`Non Proliferative`, `Proliferative`,  `Apoptotic`, `Enlarged`), names_to = 'Cell Classification', values_to = 'Percentage') 

Barplot_Average_Percentage <-ggplot(Percentage_scaled_pivot, aes(x=factor(Target, level=names), y=Percentage, fill=`Cell Classification`)) + 
  geom_bar(position="stack", stat="identity", width=.7)+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), 
                                                                axis.text.y = element_text(size=16),
                                                                axis.title.y = element_text(size=16),
                                                                legend.title = element_text(size=16),
                                                                legend.text = element_text(size=16), 
                                                                plot.title = element_text(size=16),
                                                                panel.background = element_rect(fill="white")) + 
  labs(title='', y='Average Percentage', x='')+ylim(0,101) +
  scale_fill_manual(values=c('#0072B2', '#E69F00',"#56B4E9", '#009E73'))  +
  geom_vline(xintercept = 3.5, color = 'gray30')+
  geom_vline(xintercept = 6.5, color = 'gray30')+
  geom_vline(xintercept = 9.5, color = 'gray30')+
  geom_vline(xintercept = 12.5, color = 'gray30')+
  geom_vline(xintercept = 15.5, color = 'gray30')+
  geom_vline(xintercept = 18.5, color = 'gray30')+
  geom_vline(xintercept = 21.5, color = 'gray30')+
  geom_vline(xintercept = 24.5, color = 'gray30')+
  geom_vline(xintercept = 27.5, color = 'gray30')+
  geom_vline(xintercept = 30.5, color = 'gray30')+
  geom_vline(xintercept = 33.5, color = 'gray30')


## Save file
ggsave(file.path(top_dir, 'PLOTS', 'barplot_of_average_classified_cells_scaled_data.png'), Barplot_Average_Percentage , dpi = 300, width = 12, height = 10)

#### Remove Clinical Targets - plot both ways, with and without clinical ---------------
Average_percentage_scaled_Top_Hits_only <- Average_percentage_scaled %>% filter(!Target %in% c("EAF1","EAF2","EAF1|EAF2" ,"INTS6","INTS6L","INTS6|INTS6L","TTC7A","TTC7B","TTC7A|TTC7B"))

colnames(Average_percentage_scaled_Top_Hits_only)
Percentage_scaled_Top_Hits_only<- Average_percentage_scaled_Top_Hits_only %>% group_by(Target) %>% mutate(`Non Proliferative` = mean(Percentage_Non_proliferative),
                                                                                                          `Proliferative` = mean(Percentage_Proliferative),
                                                                                                          `Apoptotic` = mean(Percentage_Apoptotic),
                                                                                                          `Enlarged` = mean(Percentage_Enlarged)) %>%
  select(c("Target", "Group_Target", "Non Proliferative" , "Proliferative", "Apoptotic", "Enlarged")) %>% unique()

Percentage_scaled_pivot_Top_Hits_only <- Percentage_scaled_Top_Hits_only %>% pivot_longer(c(`Non Proliferative`, `Proliferative`,  `Apoptotic`, `Enlarged`), names_to = 'Cell Classification', values_to = 'Percentage') 

Barplot_Average_Percentage_Top_Hits_only <-ggplot(Percentage_scaled_pivot_Top_Hits_only, aes(x=factor(Target, level=names), y=Percentage, fill=`Cell Classification`)) + 
  geom_bar(position="stack", stat="identity", width=.7)+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), 
                                                                axis.text.y = element_text(size=16),
                                                                axis.title.y = element_text(size=16),
                                                                legend.title = element_text(size=16),
                                                                legend.text = element_text(size=16), 
                                                                plot.title = element_text(size=16),
                                                                panel.background = element_rect(fill="white")) + 
  labs(title='', y='Average Percentage', x='')+ylim(0,101) +
  scale_fill_manual(values=c('#0072B2', '#E69F00',"#56B4E9", '#009E73'))  +
  geom_vline(xintercept = 3.5, color = 'gray30')+
  geom_vline(xintercept = 6.5, color = 'gray30')+
  geom_vline(xintercept = 9.5, color = 'gray30')+
  geom_vline(xintercept = 12.5, color = 'gray30')+
  geom_vline(xintercept = 15.5, color = 'gray30')+
  geom_vline(xintercept = 18.5, color = 'gray30')+
  geom_vline(xintercept = 21.5, color = 'gray30')+
  geom_vline(xintercept = 24.5, color = 'gray30')+
  geom_vline(xintercept = 27.5, color = 'gray30')+
  geom_vline(xintercept = 30.5, color = 'gray30')+
  geom_vline(xintercept = 33.5, color = 'gray30')

## Save file
ggsave(file.path(top_dir, 'PLOTS', 'barplot_of_average_classified_cells_scaled_data_top_hits_only.png'), Barplot_Average_Percentage_Top_Hits_only , dpi = 300, width = 12, height = 10)


######## Work out change in fold change for second figure - firstly plot this was all targets, then plot again removing the clinical targets-------------------------------
#####################fold change for apoptosis----------------------------------------------------------------------
#Rational - compare change in apoptosis to control only for each replicate to remove total staining batch affects
#select only apoptosis and relevant columns (note should we do this on the scaled or raw data? -checked, results exactly the same)
replicate<- scaled_objects.wide %>% select(c(`Target`, `Group_Target`, `Replicate`, `scaled_Apoptotic - Number of Objects`,  `scaled_Cells Final - Number of Objects`))

#Divide the number of apoptotic cells by the total number of cells per row
Fraction_of_apop<- replicate %>% mutate('Fraction' = (`scaled_Apoptotic - Number of Objects`/`scaled_Cells Final - Number of Objects`))

#Work out the average fraction of apoptotic cells in the controls for each replicate (view number)
N1_apop_control <- Fraction_of_apop %>% filter(Replicate == 'N1') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N1 = mean(Fraction))
N2_apop_control <- Fraction_of_apop %>% filter(Replicate == 'N2') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N2 = mean(Fraction))
N3_apop_control <- Fraction_of_apop %>% filter(Replicate == 'N3') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N3 = mean(Fraction))
N4_apop_control <- Fraction_of_apop %>% filter(Replicate == 'N4') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N4 = mean(Fraction))

print(N1_apop_control$Control_N1 %>% unique())
print(N2_apop_control$Control_N2 %>% unique())
print(N3_apop_control$Control_N3 %>% unique())
print(N4_apop_control$Control_N4 %>% unique())

N1_apop<- Fraction_of_apop %>% filter(Replicate == 'N1') %>% mutate(FoldChange_N1 =Fraction/0.234406)
N2_apop<- Fraction_of_apop %>% filter(Replicate == 'N2') %>% mutate(FoldChange_N2 =Fraction/0.1934016)
N3_apop<- Fraction_of_apop %>% filter(Replicate == 'N3') %>% mutate(FoldChange_N3 =Fraction/0.05992583)
N4_apop<- Fraction_of_apop %>% filter(Replicate == 'N4') %>% mutate(FoldChange_N4 =Fraction/0.02472507)


#Combine the fold change data for apoptosis per replicate 
FC_1_apop <- full_join(N1_apop, N2_apop)
FC_2_apop <- full_join(N3_apop, N4_apop)
Fold_Change_Apop <- full_join(FC_1_apop, FC_2_apop) %>% pivot_longer(c(FoldChange_N1, FoldChange_N2, FoldChange_N3, FoldChange_N4), names_to = 'repeat', values_to='Fold_Change') %>% drop_na() %>% filter(!Target=='Parental')


##### Repeat above code for proliferation
replicate_prolif<- scaled_objects.wide %>% select(c(`Target`, `Group_Target`, `Replicate`, `scaled_Proliferative - Number of Objects`,  `scaled_Cells Final - Number of Objects`))
Fraction_of_prolif<- replicate_prolif %>% mutate('Fraction' = (`scaled_Proliferative - Number of Objects`/`scaled_Cells Final - Number of Objects`))
N1_prol_control <- Fraction_of_prolif %>% filter(Replicate == 'N1') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N1 = mean(Fraction))
N2_prol_control <- Fraction_of_prolif %>% filter(Replicate == 'N2') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N2 = mean(Fraction))
N3_prol_control <- Fraction_of_prolif %>% filter(Replicate == 'N3') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N3 = mean(Fraction))
N4_prol_control <- Fraction_of_prolif%>% filter(Replicate == 'N4') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N4 = mean(Fraction))

print(N1_prol_control$Control_N1 %>% unique())
print(N2_prol_control$Control_N2 %>% unique())
print(N3_prol_control$Control_N3 %>% unique())
print(N4_prol_control$Control_N4 %>% unique())

N1_prolif<- Fraction_of_prolif %>% filter(Replicate == 'N1') %>% mutate(FC_N1 =Fraction/0.3317886)
N2_prolif<-Fraction_of_prolif %>% filter(Replicate == 'N2') %>% mutate(FC_N2 =Fraction/0.3511531)
N3_prolif<- Fraction_of_prolif %>% filter(Replicate == 'N3') %>% mutate(FC_N3 =Fraction/0.392028)
N4_prolif<- Fraction_of_prolif %>% filter(Replicate == 'N4') %>% mutate(FC_N4 =Fraction/0.3904623)

FC_prolif_1 <- full_join(N1_prolif, N2_prolif)
FC_prolif_2 <- full_join(N3_prolif, N4_prolif)
Fold_Change_Prolif <- full_join(FC_prolif_1, FC_prolif_2) %>% pivot_longer(c(FC_N1, FC_N2, FC_N3, FC_N4), names_to = 'repeat', values_to='Fold_Change') %>% drop_na() %>% filter(!Target=='Parental')



##### Fold change for Non-Proliferation
replicate_N_prolif<- scaled_objects.wide %>% select(c(`Target`, `Group_Target`, `Replicate`, `scaled_Non-Proliferative - Number of Objects`,  `scaled_Cells Final - Number of Objects`))
Fraction_of_N_prolif<- replicate_N_prolif %>% mutate('Fraction' = (`scaled_Non-Proliferative - Number of Objects`/`scaled_Cells Final - Number of Objects`))
N1_N_prol_control <- Fraction_of_N_prolif%>% filter(Replicate == 'N1') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N1 = mean(Fraction))
N2_N_prol_control <- Fraction_of_N_prolif %>% filter(Replicate == 'N2') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N2 = mean(Fraction))
N3_N_prol_control <- Fraction_of_N_prolif %>% filter(Replicate == 'N3') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N3 = mean(Fraction))
N4_N_prol_control <- Fraction_of_N_prolif%>% filter(Replicate == 'N4') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N4 = mean(Fraction))

print(N1_N_prol_control$Control_N1 %>% unique())
print(N2_N_prol_control$Control_N2 %>% unique())
print(N3_N_prol_control$Control_N3 %>% unique())
print(N4_N_prol_control$Control_N4 %>% unique())

N1_N_prolif<- Fraction_of_N_prolif %>% filter(Replicate == 'N1') %>% mutate(FC_N1 =Fraction/0.3961827)
N2_N_prolif<-Fraction_of_N_prolif %>% filter(Replicate == 'N2') %>% mutate(FC_N2 =Fraction/0.4166029)
N3_N_prolif<- Fraction_of_N_prolif %>% filter(Replicate == 'N3') %>% mutate(FC_N3 =Fraction/0.5142401)
N4_N_prolif<- Fraction_of_N_prolif %>% filter(Replicate == 'N4') %>% mutate(FC_N4 =Fraction/0.5550785)

FC_N_prolif_1 <- full_join(N1_N_prolif, N2_N_prolif)
FC_N_prolif_2 <- full_join(N3_N_prolif, N4_N_prolif)
Fold_Change_N_Prolif <- full_join(FC_N_prolif_1, FC_N_prolif_2) %>% pivot_longer(c(FC_N1, FC_N2, FC_N3, FC_N4), names_to = 'repeat', values_to='Fold_Change') %>% drop_na()%>% filter(!Target=='Parental')


##### Fold change for Enlarged
replicate_Enlarged<- scaled_objects.wide %>% select(c(`Target`, `Group_Target`, `Replicate`, `scaled_Enlarged - Number of Objects`,  `scaled_Cells Final - Number of Objects`))
Fraction_of_Enlarged<- replicate_Enlarged %>% mutate('Fraction' = (`scaled_Enlarged - Number of Objects`/`scaled_Cells Final - Number of Objects`))
N1_Enlarged_control <- Fraction_of_Enlarged%>% filter(Replicate == 'N1') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N1 = mean(Fraction))
N2_Enlarged_control <- Fraction_of_Enlarged %>% filter(Replicate == 'N2') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N2 = mean(Fraction))
N3_Enlarged_control <- Fraction_of_Enlarged %>% filter(Replicate == 'N3') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N3 = mean(Fraction))
N4_Enlarged_control <- Fraction_of_Enlarged%>% filter(Replicate == 'N4') %>% filter(Target %in% c('Control_1', 'Control_2', 'Control_1|Control_2')) %>% mutate(Control_N4 = mean(Fraction))

print(N1_Enlarged_control$Control_N1 %>% unique())
print(N2_Enlarged_control$Control_N2 %>% unique())
print(N3_Enlarged_control$Control_N3 %>% unique())
print(N4_Enlarged_control$Control_N4 %>% unique())

N1_Enlarged<- Fraction_of_Enlarged %>% filter(Replicate == 'N1') %>% mutate(FC_N1 =Fraction/0.002673601)
N2_Enlarged<-Fraction_of_Enlarged %>% filter(Replicate == 'N2') %>% mutate(FC_N2 =Fraction/0.003971227)
N3_Enlarged<- Fraction_of_Enlarged %>% filter(Replicate == 'N3') %>% mutate(FC_N3 =Fraction/0.003382909)
N4_Enlarged<- Fraction_of_Enlarged %>% filter(Replicate == 'N4') %>% mutate(FC_N4 =Fraction/0.008309519)

FC_Enlarged_1 <- full_join(N1_Enlarged, N2_Enlarged)
FC_Enlarged_2 <- full_join(N3_Enlarged, N4_Enlarged)
Fold_Change_Enlarged <- full_join(FC_Enlarged_1, FC_Enlarged_2) %>% pivot_longer(c(FC_N1, FC_N2, FC_N3, FC_N4), names_to = 'repeat', values_to='Fold_Change') %>% drop_na() %>% filter(!Target=='Parental')

####Figure format for paper
Fold_change_apop_combined<- ggplot(Fold_Change_Apop, aes(x=factor(Target, level=names), y=Fold_Change, fill=Group_Target)) + 
  geom_boxplot() + theme_classic()+  theme(axis.text.x = element_blank(), 
                                           legend.position = 'none',
                                           axis.text.y = element_text(size=16),
                                           axis.title.y = element_text(hjust=0.5,size=16),
                                           plot.title = element_text(size=16)) + 
  labs(title='Apoptotic', y='Fold Change in \nApoptotic Cells', x='') 

Fold_change_Enlarged_combined <- ggplot(Fold_Change_Enlarged, aes(x=factor(Target, level=names), y=Fold_Change, fill=Group_Target)) + 
  geom_boxplot() +  theme_classic()+  theme(axis.text.x = element_blank(), 
                                            legend.position = 'none',
                                            axis.text.y = element_text(size=16),
                                            axis.title.y = element_text(hjust=0.5, size=16),
                                            plot.title = element_text(size=16)) +
  labs(title='Enlarged', y='Fold Change in \nEnlarged Cells', x='')


Fold_change_N_prolif_combined <- ggplot(Fold_Change_N_Prolif, aes(x=factor(Target, level=names), y=Fold_Change, fill=Group_Target)) + 
  geom_boxplot() +  theme_classic()+  theme(axis.text.x = element_blank(), 
                                            legend.position = 'none',
                                            axis.text.y = element_text(size=16),
                                            axis.title.y = element_text(hjust=0.5,size=16),
                                            plot.title = element_text(size=16)) + 
  labs(title='Non Proliferative', y='Fold Change in \nNon-proliferative Cells', x='')



Fold_change_prolif_combined <- ggplot(Fold_Change_Prolif, aes(x=factor(Target, level=names), y=Fold_Change, fill=Group_Target)) + 
  geom_boxplot() +  theme_classic()+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16),
                                            legend.position = 'none',
                                            axis.text.y = element_text(size=16),
                                            axis.title.y = element_text(hjust=0.5,size=16),
                                            plot.title = element_text(size=16)) + 
  labs(title='Proliferative', y='Fold Change in \nProliferative Cells', x='')

#Arrange scaled fold change compared to control average plots
Fold_Change_Classified_Cells <- grid.arrange(Fold_change_apop_combined,
                                             Fold_change_Enlarged_combined,
                                             Fold_change_N_prolif_combined,
                                             Fold_change_prolif_combined,
                                             nrow=4,
                                             heights = c(1, 1, 1, 1.5))

## Save file
ggsave(file.path(top_dir, 'PLOTS', 'Fold_Change_Classified_Cells_Scaled_FC_to_Averaged_controls_minus_parentals.png'), 
       Fold_Change_Classified_Cells , dpi = 300, width = 12, height = 18)

#### This figure still needs the stats adding