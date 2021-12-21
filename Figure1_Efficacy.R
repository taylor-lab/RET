require(data.table)
require(PieDonut)
require(ggplot2)
require(RColorBrewer)
require(wesanderson)
require(grid)
require(gridExtra)
require(moonBook)
require(webr)


source('data/grob_stuff.R')

data <- fread('data/RET_master_table.tsv', header = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 1a/d: Pie Chart & Overall Efficacy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##### Pie Chart Figure 1a #####

PieDonut(data,aes(pies=Mut_Type,donuts=Tumor_Type), labelposition=2, explode = 1, labelpositionThreshold=0.5)

##### Waterfall Patient Order ######

recist_order = unique(data[,.(`Study_ID`,`RECIST_response`)])[order(`RECIST_response`, decreasing=T)]$`Study_ID`

##### Efficacy Figure 1d #####

do.waterfall = function(br, recist_order, leg.pos='none'){
  
  w.f = ggplot(br,aes(x=factor(`Study_ID`, recist_order),
                      y=`RECIST_response`)) +
    geom_bar(stat='identity', width=.7, fill = 'black') +
    labs(x='',y='% best Change') +
    theme_classic(base_size=10) +
    ylim(-115,110) + 
    geom_point(data=br[`Target Lesion Response`=='NTL'], aes(x=factor(`Study_ID`, recist_order),y=0),shape=5,size=2.5 ,colour='black') +
    geom_point(data=br[`Target Lesion Response`=='No Assessments'], aes(x=factor(`Study_ID`, recist_order),y=0),shape=18,size=4.5 ,colour='black') +
    geom_hline(yintercept = 0, color = 'black', size = .4) +
    geom_hline(yintercept = -30, color = '#DB4512', size = .5, linetype='dashed') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.y = element_line(color="black", size=.6),
          legend.key.size=unit(0.35,'cm'),
          legend.position=leg.pos,
          legend.key = element_blank(),
          legend.background = element_blank(),
          axis.line.x=element_blank()) +
    guides(colour = guide_legend(override.aes = list(shape = NA)))
  
  w.f
  
}

w.f = do.waterfall(data, recist_order)

#Alteration Type#
do.alttype = function(at, recist_order, leg.pos='none'){
  
  a.t = ggplot(at) + geom_point(aes(x=factor(Study_ID, recist_order),y=1, colour=`Alt_Bin_Det`), shape=15, size=3) + 
    scale_colour_manual(values=c('KD_oncogenic'='#9CB770','CRD_oncogenic'='#FCB89B',
                                 'CRD_indel_oncogenic'='#FDD29B','KD_pathogenic'='#69A583',
                                 'CRD_pathogenic'='#EA7E82', 'TMD VUS'='grey',
                                 'CCDC6'='#FF4E50', 'KIF5B'='#FFCA06', 'NCOA4'='#8C4B21',
                                 'TAF3'='#B49A85', 'RUFY3'='#B49A85', 'KIF13A'='#B49A85', 'ERC1'='#B49A85',
                                 'EML4'='#B49A85', 'KIAA1468'='#B49A85', 'CLIP1'='#B49A85', 'Unknown'='#8A9294',
                                 'RBPMS'='#B49A85'))  +
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
  
  a.t
}

a.t = do.alttype(data, recist_order)

#Tumor Type#
do.tumortype = function(tt, recist_order, leg.pos='none'){
  
  t.t = ggplot(tt) + geom_point(aes(x=factor(Study_ID, recist_order),y=1, colour=`Tumor_Type`), shape=15, size=3) + 
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
  
  t.t
}

t.t = do.tumortype(data, recist_order)

#Overall Response#
do.overallchange = function(oc, recist_order, leg.pos='none'){
  
  my.cols = brewer.pal(8,'Dark2')
  
  o.c = ggplot(oc) + geom_point(aes(x=factor(Study_ID, recist_order),y=1, colour=`Overall Response`), shape=15, size=3) + 
    scale_colour_manual(values=c('CR'='firebrick2', 'PR'='indianred1', 'SD'='#278c44', 'PD'='#0072B2','NE'='grey'),
                        breaks=c('CR'='firebrick2', 'PR'='indianred1', 'SD'='#278c44', 'PD'='#0072B2','NE'='grey'),
                        labels=c('CR'='firebrick2', 'PR'='indianred1', 'SD'='#278c44', 'PD'='#0072B2','NE'='grey'),
                        name='') +
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
  
  o.c
}

o.c = do.overallchange(data, recist_order)

#Mutation Origin#
do.conc = function(conc, recist_order, leg.pos='none'){
  
  c.f = ggplot(conc) + geom_point(aes(x=factor(Study_ID, recist_order),y=1,colour=Mut_Type),shape=15, size=3) +
    scale_colour_manual(values=c('Germline Mutation'='black','Somatic Mutation'='grey', 'Fusion'='grey')) +
    scale_y_continuous(limits=c(1,1)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(angle=0, hjust=.5, size=10),
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
  
  c.f
}

c.f = do.conc(data, recist_order)

## Previous RET TKI ##
do.prev = function(prev, recist_order, leg.pos='none'){
  
  p.v = ggplot(prev) + geom_point(aes(x=factor(Study_ID, recist_order),y=1,colour=Previous_TKI),shape=15, size=3) +
    scale_colour_manual(values=c('Yes'='black','No'='grey')) +
    scale_y_continuous(limits=c(1,1)) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(angle=0, hjust=.5, size=10),
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
  
  p.v
}

p.v = do.prev(data, recist_order)

#Duration#
do.duration = function(dr, recist_order, leg.pos='none'){
  
  my.cols = brewer.pal(9,'Blues')
  
  o.d = ggplot(dr,aes(x=factor(Study_ID, recist_order),
                      y=Weeks)) + 
    geom_bar(stat='identity',width=.7, fill='grey') +
    labs(x='',y='Duration (wks)') +
    theme_classic(base_size=10) +
    geom_hline(yintercept = 0, color = 'black', size = .4) +
    ylim(0,130) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.y = element_line(color="black", size=.7),
          legend.key.size=unit(0.35,'cm'),
          legend.position=leg.pos,
          legend.key = element_blank(),
          legend.background = element_blank(),
          axis.line.x=element_blank()) +
    geom_point(data=dr[`Ongoing`=='YES'],aes(x=factor(Study_ID, recist_order),y=Weeks+1),shape=17,size=2.4 ,colour='grey') +
    guides(colour = guide_legend(override.aes = list(shape = NA)))
  
  o.d
  
}

o.d = do.duration(data, recist_order)

w.f.f = w.f + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
a.t.f = a.t + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
t.t.f = t.t + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
o.c.f = o.c + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
c.f.f = c.f + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
p.v.f = p.v + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())
o.d.f = o.d + facet_grid(.~Alt_Type,scales='free_x',space='free_x') + theme(strip.background = element_blank(), strip.text = element_blank())

p.grobs = align.grobs(list(w.f.f, a.t.f, t.t.f, o.c.f,  c.f.f, p.v.f, o.d.f))
c.plot <- grid.arrange(arrangeGrob(p.grobs[[1]],p.grobs[[2]],p.grobs[[3]],p.grobs[[4]],p.grobs[[5]],p.grobs[[6]],p.grobs[[7]],
                                   heights=c(.5, .1, .1, .1, .1, .1, .2)))
