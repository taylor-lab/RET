require(data.table)
require(PieDonut)
require(ggplot2)
require(grid)
require(gridExtra)
require(dplyr)
require(Rmisc)

source('data/grob_stuff.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 4a/b/c: Whole Exome Sequencing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### Data #####
data <- fread('data/RET_master_table.tsv', header = T)
exome = fread('data/EDfig4a_data.tsv', header = T)
exome_ids = fread('data/exome_cases.tsv', header = F)
tmb = fread('data/EDfig4b.tsv', header = T)

# Set Fig 4a patient order #
exome_order0 = data %>%
  filter(Study_ID %in% exome_ids$V1) %>%
  select(Study_ID, Tumor_Bin, Weeks, Ongoing, Benefit_OR, Benefit_TL) %>%
  arrange(Tumor_Bin, Benefit_TL, Weeks)

exome_order = as.character(exome_order0$Study_ID)

##### ED Fig 4a #####
tmb_plot = ggplot(exome, aes(x=factor(Study_ID, exome_order), y=TMB)) +
  geom_bar(stat='identity', width=0.7, fill='black') +
  geom_hline(yintercept = 0, color = 'black', size = .4) +
  geom_hline(yintercept = 2.9, color = '#74a089', size = .5, linetype='dashed') +
  geom_hline(yintercept = 1.03, color = '#1d4e89', size = .5, linetype='dashed') +
  geom_hline(yintercept = 1.6, color = '#9dc6d8', size = .5, linetype='dashed') +
  ylim(0,10) +
  theme_bw() +
  ylab('TMB (Mut/Mb)') +
  xlab('') +
  theme(legend.position = "none")

tmb.f = tmb.f + facet_grid(.~Type, scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
tmb.f

#Tumor Type#
do.tumortype = function(tt, exome_order, leg.pos='none'){
  
  t.t = ggplot(tt) + geom_point(aes(x=factor(Study_ID, exome_order),y=1, colour=`Tumor`), shape=15, size=3) + 
    scale_colour_manual(values=c('Papillary Thyroid'='#9dc6d8','Medullary Thyroid'='#1d4e89',
                                 'Poorly Differentiated Thyroid'='#00b3ca',
                                 'Lung Adenocarcinoma'='#74a089','Lung Adenosquamous'='#f69256',
                                 'Combined Small Cell Lung Carcinoma'='#e38690',
                                 'High-Grade Neuroendocrine Carcinoma of the Colon and Rectum'='#513c1f', 
                                 'Non Langerhans Cell Histiocytosis'='#4e5c5f')) +
    scale_y_continuous(limits=c(1,1)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          legend.key.size=unit(0.35,'cm'),
          legend.position=leg.pos,
          plot.margin = unit(c(0,0,0,0), 'lines'),
          legend.key=element_blank(),
          legend.background = element_blank())
  
  t.t
}

t.t = do.tumortype(exome, exome_order)
t.t.f = t.t + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
t.t.f


#Overall Response#
do.overallchange = function(oc, exome_order, leg.pos='none'){
  
  my.cols = brewer.pal(8,'Dark2')
  
  o.c = ggplot(oc) + geom_point(aes(x=factor(Study_ID, exome_order),y=1, colour=`Overall Response`), shape=15, size=3) + 
    scale_colour_manual(values=c('CR'='firebrick2', 'PR'='indianred1', 'SD'='#278c44')) +
    scale_y_continuous(limits=c(1,1)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          #axis.title.y=element_text(angle=0, hjust=.5, size=10), #####
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          legend.key.size=unit(0.35,'cm'),
          legend.position=leg.pos,
          plot.margin = unit(c(0,0,0,0), 'lines'),
          legend.key=element_blank(),
          legend.background = element_blank())
  
  o.c
}

o.c = do.overallchange(exome, exome_order)
o.c.f = o.c + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
o.c.f

#SMOKING TRACK#
do.smoking = function(sm, exome_order, leg.pos='none'){
  
  s.m = ggplot(sm) + geom_point(aes(x=factor(Study_ID, exome_order),y=1, colour=`Smoking`), shape=15, size=3) + 
    scale_colour_manual(values=c('yes'='black',
                                 'no'='white',
                                 'NE'='white')) +
    scale_y_continuous(limits=c(1,1)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          #axis.title.y=element_text(angle=0, hjust=.5, size=10), #####
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          legend.key.size=unit(0.35,'cm'),
          legend.position=leg.pos,
          plot.margin = unit(c(0,0,0,0), 'lines'),
          legend.key=element_blank(),
          legend.background = element_blank())
  
  s.m
}

s.m = do.smoking(exome, exome_order)
s.m.f = s.m + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
s.m.f


sig0 = tbl_df(exome) %>% select(-Tumor, -Alt_Type, -Ongoing,  -Smoking, -`Overall Response`, -TMB) %>%
  reshape2::melt(id.vars = 'Study_ID') %>%
  mutate(variable = stringr::str_replace(variable, 'Signature_', ''))

sig0 = mutate(sig0, variable = mapvalues(variable, seq(1,31),
                                         c('Age','APOBEC','BRCA','Smoking','Other','MMR','Other','Other','Other','Other','Other','Other','APOBEC','Other','MMR',
                                           'Other','Other','Other','Other','MMR','Other','Other','Other','Other','Other','MMR','Other','Other','Other','Other', 'Other'))) %>%
  mutate(variable = factor(variable, levels = c('Age','APOBEC','MMR','BRCA','Smoking','Other'), ordered = T))

col_vec = c("#939598","#FBB040","#BE1E2D", "#AE742A", "#FFDCAE", "#000000")

sign_plot = ggplot(sig0, aes(x = factor(`Study_ID`, exome_order), fill = value)) +
  geom_bar(aes(weight = value, fill = variable), width=.7) +
  theme_bw() +
  scale_fill_manual(values = col_vec) +
  theme(text = element_text(size = 12), legend.key = element_blank()) +
  guides(fill = guide_legend(keywidth = .5, keyheight = .5, ncol = 2, title = 'Signature')) +
  ylab('Fraction of mutations') +
  xlab('') +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "bottom")

sign_plot = sign_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#(labels == FALSE) sign_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


sign.f = sign_plot + facet_grid(.~exome$Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
sign.f

p.grobs = align.grobs(list(tmb.f, t.t.f, o.c.f, s.m.f, sign.f))
c.plot <- grid.arrange(arrangeGrob(p.grobs[[1]],p.grobs[[2]],p.grobs[[3]],p.grobs[[4]],p.grobs[[5]],
                                   heights=c(.3, .1, .1, .1, .4)))


###### ED Fig 4b #####

tmb = tmb %>%
  filter(RECIST!='NE')
tmb$quartile <- ntile(tmb$TMB, 4)
tmb$RECIST <- as.numeric(tmb$RECIST)

tmb1 = summarySE(tmb, measurevar="TMB", groupvars=c("quartile"))

box3 = ggplot(tmb1, aes(x = as.character(quartile), y = TMB)) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, color="red") +
  geom_errorbar(aes(ymin=TMB-se, ymax=TMB+se), colour="black", width=.1) +
  ylab('TMB (Mut/Mb)') +
  xlab('Percentile of TMB') +
  theme_bw()

box3

box2 = ggplot(tmb, aes(x = as.character(quartile), y = RECIST)) +
  geom_boxplot(alpha = 0.60) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  geom_hline(yintercept = -30, color = '#DB4512', size = .5, linetype='solid') +
  ylab('% best change') +
  xlab('Percentile of TMB') +
  ylim(-100,0) +
  theme_bw()

box2

anova_one_way <- aov(RECIST~quartile, data = tmb)
summary(anova_one_way) #0.422

box1 = ggplot(tmb, aes(x = as.character(quartile), y = Weeks)) +
  geom_boxplot(alpha = 0.60) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  ylab('Weeks on therapy') +
  xlab('Percentile of TMB') +
  theme_bw()
box1

anova_one_way <- aov(Weeks~quartile, data = tmb)
summary(anova_one_way) #0.479

p.grobs = align.grobs(list(box3, box2, box1))
c.plot <- grid.arrange(arrangeGrob(p.grobs[[1]],p.grobs[[2]],p.grobs[[3]],
                                   heights=c(.5, .5, .5)))


##### ED Fig 4c #####
box5 = ggplot(tmb, aes(x = Clonality, y = Weeks)) +
  geom_boxplot(alpha = 0.60) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  ylab('Weeks on therapy') +
  xlab('Degree of subclonality') +
  theme_bw()
box5

anova_one_way <- aov(Weeks~Clonality, data = tmb)
summary(anova_one_way) #0.976

box6 = ggplot(tmb, aes(x = Clonality, y = RECIST)) +
  geom_boxplot(alpha = 0.60) +
  geom_jitter(color="grey", size=2, alpha=0.9) +
  geom_hline(yintercept = -30, color = '#DB4512', size = .5, linetype='solid') +
  ylab('% best change') +
  xlab('Degree of subclonality') +
  ylim(-100,0) +
  theme_bw()
box6

anova_one_way <- aov(RECIST~Clonality, data = tmb)
summary(anova_one_way) #0.68

p.grobs = align.grobs(list(box6, box5))
c.plot <- grid.arrange(arrangeGrob(p.grobs[[1]],p.grobs[[2]],
                                   heights=c(.5, .5)))