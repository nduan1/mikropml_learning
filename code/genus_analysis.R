library(tidyverse)
library(broom)
library(ggtext)
shared <- read_tsv("raw_data/baxter.subsample.shared",
         col_types = cols(Group=col_character(),
                          .default = col_double())) %>% 
  rename_all(tolower) %>% #make all col name lower case
  select(group,starts_with("otu")) %>% 
  pivot_longer(-group,names_to = "otu",values_to = "count")

taxonomy <- read_tsv("raw_data/baxter.cons.taxonomy") %>% 
  rename_all(tolower) %>% 
  select(otu,taxonomy) %>% #select cols
  mutate(otu=tolower(otu),
         taxonomy=str_replace_all(taxonomy,"\\(\\d+\\)",""),
         taxonomy=str_replace(taxonomy,";unclassified","_unclassified"),#to replace the first match
         taxonomy=str_replace_all(taxonomy,";unclassified",""),
         taxonomy=str_replace_all(taxonomy,";$",""),#remove final semi colon
         taxonomy=str_replace_all(taxonomy,".*;","")#remove upper stream info before the last item 
         )  #we kept the genus name

metadata <- read_tsv("raw_data/baxter.metadata.tsv",
         col_types = cols(sample=col_character())) %>% 
  rename_all(tolower) %>% 
  rename(group=sample) %>% 
  mutate(srn=dx_bin=="Adv Adenoma" | dx_bin=="Cancer",
         lesion=dx_bin=="Adv Adenoma" | dx_bin=="Cancer"| dx_bin=="Adenoma")

###join file##
composite <- inner_join(shared,taxonomy,by="otu") %>% 
  group_by(group,taxonomy) %>% 
  summarise(count=sum(count), .groups = "drop") %>% 
  group_by(group) %>% #group is patient id
  mutate(rel_abund=count/sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  inner_join(.,metadata, by="group")

###cc122 Figure out which genera are significnatly different between with or without disease
sig_genera <- composite %>% 
  nest(data=-taxonomy) %>% #data is the columns wihtout taxonomy
  mutate(test=map(.x=data,~ wilcox.test(rel_abund~srn,data =.x ) %>% tidy)) %>% #map help iterate over different rows
  unnest(test) %>% 
  mutate(p.adj=p.adjust(p.value,method = "BH")) %>% 
  filter(p.adj<0.05) %>% 
  select(taxonomy,p.adj)

composite %>% 
  inner_join(sig_genera,by="taxonomy") %>% 
  mutate(rel_abund=100* (rel_abund+1/2000),
         taxonomy=str_replace(taxonomy,"(.*)","*\\1*"),
         taxonomy=str_replace(taxonomy,"\\*(.*)_unclassified\\*","Unclassified<br>*\\1*"),
         srn=factor(srn,levels =c(T,F) )) %>% #\\1 means take everything in ()
  ggplot(aes(x=rel_abund,y=taxonomy,color=srn))+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.3),
              shape=21)+
  stat_summary(fun.data = median_hilow,fun.args = list(conf.int=0.5),
               geom="pointrange",
               position = position_dodge(width=0.8),
               color="black",show.legend = F)+
  scale_x_log10()+
  scale_color_manual(NULL,
                     breaks = c(F, T),
                     values=c("gray","dodgerblue"),#set up legend order
                     labels=c("Healthy","SRN"))+
  labs("Relative abundance(%)",y=NULL)+
  theme_classic()+
  theme(
    axis.text.y=element_markdown()
  )

ggsave("figures/significant_genera.tiff",width = 6,height = 4)
###statistical test not tell you if the genus associate with certain 
#the 0 median value of all these genera indicate that more than half of the 
#people don't have these bacteria. there could be a bacterial population that more poorly 
#distributed (e.g. like within 10% people) or could be a stronger biomarkers of srn
#if 10% people have that marker, that is the sure sign they can get srn. But statistical
#test won't tell you which one is biomarker.








