
df = read.csv("F:/ZLF/CAP/data/relative_data/dfall-220617.csv",
              row.names=1)
metadata = read.csv('F:/ZLF/CAP/data/relative_data/metadata-all-220617.csv',row.names = 1)



# PERMANOVA ---------------------------------------------------------------
vals = c('sex','age','BMI', 'city',
         'smoke','basic_LRT_disease','pre_biotic','pre_immunosuppressive',
         'respiratory.support.invasive','death_case1','severe_case') 

sub.meta = filter(metadata,d == 1)
for (val in test.vals) {sub.meta <- sub.meta[!is.na(sub.meta[,val]),]}
sub.df = df[,rownames(sub.meta)] %>% t() %>% as.data.frame()

OTU = phyloseq::otu_table(sub.df, taxa_are_rows = F)
physeq = phyloseq(OTU)
jsd=phyloseq::distance(physeq, method = "jsd")
res <- adonis2(as.formula(paste0("jsd~",'sex+age+BMI+city+basic_LRT_disease+respiratory.support.invasive+pre_immunosuppressive+pre_biotic+smoke+death_case1+severe_case')),data = sub.meta,permutations = 999,by = 'margin')
res$'adjust.p' <- p.adjust (res$`Pr(>F)`, method = 'fdr', n = (length(vals)))

res


multi_res = data.frame(R2 = rep(NA,length(vals)),
                       pval = rep(NA,length(vals))
)
row.names(multi_res) <- vals

multi_res$R2 <- res[vals,'R2']
multi_res$pval <- res[vals,'Pr(>F)']
multi_res$padj <- res[vals,'adjust.p']

multi_res$group <- 'multi'
multi_res <- multi_res[order(multi_res[,'R2'],decreasing = T),]
multi_res$val <- rownames(multi_res)

multi_res

# cluster -----------------------------------------------------------------
## PAM
pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))

max_num = 20

nums=c(2:{max_num})
cord=df %>% t() %>% as.data.frame()

for(num in nums){{
  pam = pam1(cord,k=num)
  ofile=c("cord_pam_",num,".csv")
  write.csv(pam ,ofile)    
}}


# data table-----------------------------------------------------------------

vals = c('sex','city',
         'smoke','basic_LRT_disease','pre_biotic','pre_immunosuppressive',
         'respiratory.support.invasive','death_case1','severe_case') 

vals = c('age','BMI') 
m1 <- metadata %>%filter(subject != 'nc') %>% filter(subject != 'NCPCR') %>% filter(subject != 'healthy') %>%  t() %>% as.data.frame() %>% arrange() %>% t() %>% as.data.frame() 
m1$dfirst1 <- NA
m1$dfirst1[1] = 'dfirst'
subj <- group_by(m1,subject) %>% summarise(count = n()) %>% as.data.frame() %>% select(count)
subj1 <- subj$count

su1 = 1
for (su in subj1) {
  m1$dfirst1[su1] <- 'dfirst' 
  su1 = su1+su
}

vals = c('sex','city',
         'smoke','basic_LRT_disease','pre_biotic','pre_immunosuppressive',
         'respiratory.support.invasive','death_case1','severe_case') 

sub.meta = m1
for (val in vals) {
  indi <- sub.meta[!is.na(sub.meta[,val]),val] %>% unique()
  assign(paste0(val,'1'),sub.meta[!is.na(sub.meta[,val]),] %>% select(subject) %>% unique() %>% nrow())
  print(paste0('######',val,'######'))
  for (ind in indi) {
    assign(paste0(val,'_',ind),sub.meta[!is.na(sub.meta[,val]),] %>% filter(get(val) == ind) %>% select(subject) %>% unique() %>% nrow())
    print(paste0(paste0(val,'_',ind),'=',get(paste0(val,'_',ind)),'  ',val,'=',get(paste0(val,'1')),'  ','precentage=',get(paste0(val,'_',ind))/get(paste0(val,'1'))))
  }
  
}

sub.meta = m1 %>% filter(d == 1)
for (val in vals) {
  indi <- sub.meta[!is.na(sub.meta[,val]),val] %>% unique()
  assign(paste0(val,'1'),sub.meta[!is.na(sub.meta[,val]),] %>% select(subject) %>% unique() %>% nrow())
  print(paste0('######',val,'######'))
  for (ind in indi) {
    assign(paste0(val,'_',ind),sub.meta[!is.na(sub.meta[,val]),] %>% filter(get(val) == ind) %>% select(subject) %>% unique() %>% nrow())
    print(paste0(paste0(val,'_',ind),'=',get(paste0(val,'_',ind)),'  ',val,'=',get(paste0(val,'1')),'  ','precentage=',get(paste0(val,'_',ind))/get(paste0(val,'1'))))
  }
  
}

vals = c('age','BMI') 
for (val in vals) {
  print(val)
  print(sub.meta[!is.na(sub.meta[,val]),] %>% select(subject) %>% unique() %>% nrow())
  print(sub.meta[!is.na(sub.meta[,val]),val] %>% as.numeric() %>% quantile())
}


i
for (i in 1:10) {
  a <- filter(metadata,pam_10_cluster == i) %>% nrow()
  print(paste0(i,'  ',a))
}



# CS transfer test ----------------------------------------------------------------

transmatrix = read.csv('F:/ZLF/CAP/paper_structure/figure3/markovchain_result/ts_trans_counting_matrix_wpl.csv',row.names = 1)
transmatrix = Cluster transmatrix
#View(transmatrix)
colnames(transmatrix)<-c('1','10','2','3','4','5','6','7','8','9')
rownames(transmatrix)<-c('1','10','2','3','4','5','6','7','8','9')

clus<-c('1','10','2','3','4','5','6','7','8','9')


pva<-data.frame(cs = 1,trans_to_cs = 1,pvalue = NA)

for (k in as.character(1:10)) {
  for (i in as.character(1:10)) {
    #assign(paste0('trans',j,i))
    trans<-data.frame(matrix(nrow = 2,ncol = 2))
    trans<-data.frame(matrix(nrow = 2,ncol = 2))
    rownames(trans)<-c('cs','csothers')
    colnames(trans)<-c('transtocs','transtoothers')
    
    
    trans[1,1]<-transmatrix[i,k]
    mid = 0
    for (j in setdiff(setdiff(clus,k),i)){mid = mid + transmatrix[i,j]}
    trans[2,1]<-mid
    mid = 0
    for (j in setdiff(setdiff(clus,k),i)){mid = mid + transmatrix[j,k]}
    trans[1,2]<-mid
    mid = 0
    for (j in setdiff(clus,i)){for (h in setdiff(setdiff(clus,k),j)) {mid = mid + transmatrix[j,h]}}
    trans[2,2]<-mid
    pva[10*(as.numeric(k)-1)+as.numeric(i),'cs'] = as.numeric(i)
    pva[10*(as.numeric(k)-1)+as.numeric(i),'trans_to_cs'] = as.numeric(k)
    pva[10*(as.numeric(k)-1)+as.numeric(i),'pvalue'] = fisher.test(trans)$p
    
  }                
}

pva$qvalue<-p.adjust(pva$pvalue,method = 'BH')
for (i in nrow(pva):1) {
  if(pva[i,'cs'] == pva[i,'trans_to_cs']){
    pva = pva[-i,]
  }
  
}
pva_sig<-filter(pva,qvalue <= 0.05)

pva_sig
#pva_sig$OR<-rep(NA,nrow(pva_sig))



# Intubation related CS in severe case ------------------------------------


transmeta = read.xlsx('F:/ZLF/CAP/data/relative_data/trans_meta1.xlsx',1,rowNames = T)
transmeta<-filter(transmeta,severity == 'yes')
transmeta_mv<-filter(transmeta,respiratory_support == '1' | respiratory_support == '2,1')
transmeta_oth<-filter(transmeta,respiratory_support != '1' & respiratory_support != '2,1')


trans<-data.frame(matrix(nrow = 2,ncol = 2))
trans<-data.frame(matrix(nrow = 2,ncol = 2))
rownames(trans)<-c('mv','others')
colnames(trans)<-c('self','to others')

trans[1,1]<-nrow(filter(transmeta_mv, d1 ==d5))
trans[1,2]<-nrow(filter(transmeta_mv, d1 !=d5))
trans[2,1]<-nrow(filter(transmeta_oth, d1 ==d5))
trans[2,2]<-nrow(filter(transmeta_oth, d1 !=d5))
fisher.test(trans)

i = 5
for (i in 1:10) {
  trans<-data.frame(matrix(nrow = 2,ncol = 2))
  trans<-data.frame(matrix(nrow = 2,ncol = 2))
  colnames(trans)<-c('cs','csothers')
  rownames(trans)<-c('mv','non-mv')
  
  trans[1,1]<-nrow(filter(transmeta_mv, d1 != i & d5 == i))
  trans[1,2]<-nrow(filter(transmeta_mv, d1 !=i & d5 != i & d5 != d1))
  trans[2,1]<-nrow(filter(transmeta_oth, d1 != i & d5 == i))
  trans[2,2]<-nrow(filter(transmeta_oth, d1 !=i & d5 != i & d5 != d1))
  print(paste('cluter',i,'pvalue = ',fisher.test(trans)$p))
}




# self transfer rate ------------------------------------------------------

transmatrix = read.csv('F:/ZLF/CAP/paper_structure/figure3/markovchain_result/n_ts_trans_counting_matrix_wpl.csv',row.names = 1)

#View(transmatrix)
colnames(transmatrix)<-c('1','10','2','3','4','5','6','7','8','9')
rownames(transmatrix)<-c('1','10','2','3','4','5','6','7','8','9')
transmatrix<-select(transmatrix,'1','2','3','4','5','6','7','8','9','10')
transmatrix<-transmatrix[c('1','2','3','4','5','6','7','8','9','10'),]

transmatrix<-select(transmatrix,'1','2','3','4','5','6','7','8','9','10')
transmatrix<-transmatrix[c('1','2','3','4','5','6','7','8','9','10'),]

clus<-c('1','10','2','3','4','5','6','7','8','9')



pva<-data.frame(cs = 1,pvalue = NA)

for (i in 1:10) {
  trans<-data.frame(matrix(nrow = 2,ncol = 2))
  trans<-data.frame(matrix(nrow = 2,ncol = 2))
  colnames(trans)<-c('self','others')
  rownames(trans)<-c('csi','csothers')
  
  trans[1,1]<-transmatrix[i,i]
  mid = 0
  for (j in setdiff(clus,i)){mid = mid + transmatrix[j,j]}
  trans[2,1]<-mid
  mid = 0
  for (j in setdiff(clus,i)){mid = mid + transmatrix[i,j]}
  trans[1,2]<-mid
  mid = 0
  for (j in setdiff(clus,i)){for (h in setdiff(setdiff(clus,i),j)) {mid = mid + transmatrix[j,h]}}
  trans[2,2]<-mid
  pva[i,'cs'] = as.numeric(i)
  pva[i,'pvalue'] = fisher.test(trans)$p
  print(paste('cluter',i,'pvalue = ',fisher.test(trans)$p))
}




pva$qvalue<-p.adjust(pva$pvalue,method = 'BH')

pva_sig<-filter(pva,qvalue <= 0.05)

pva_sig




####################DAY1-5每一个转移 

# day1-5 transfer preference ----------------------------------------------

unique(transmeta$respiratory_support)
pva<-data.frame(cs = 1,trans_to_cs = 1,pvalue = NA)
transmeta = read.xlsx('F:/ZLF/CAP/data/relative_data/trans_meta1.xlsx',1,rowNames = T)
transmeta<-filter(transmeta,severity == 'yes')
#transmeta<-filter(transmeta,respiratory_support == '1' | respiratory_support == '2,1')
transmeta<-transmeta[!is.na(transmeta[,'d1']),]
transmeta<-transmeta[!is.na(transmeta[,'d5']),]
transmeta<-filter(transmeta,severity == 'yes'|severity == 'no')
j = 8
i = 5

#View(transmeta)

#transmeta$severity

#transmeta$severity
for (j in 1:10) {
  for (i in 1:10) {
    #assign(paste0('trans',j,i))
    trans<-data.frame(matrix(nrow = 2,ncol = 2))
    trans<-data.frame(matrix(nrow = 2,ncol = 2))
    rownames(trans)<-c('cs','csothers')
    colnames(trans)<-c('transtocs','transtoothers')
    
    trans[1,1]<-nrow(filter(transmeta[!is.na(transmeta[,'d1']),],d1 == j & d5 ==i))
    trans[2,1]<-nrow(filter(transmeta[!is.na(transmeta[,'d1']),],d1 != j & d1 != i & d5 == i))
    trans[1,2]<-nrow(filter(transmeta[!is.na(transmeta[,'d1']),],d1 == j & d5 != i &d5 != j))
    mid = 0
    #for (k in setdiff(setdiff(1:10,j),i)) {for(h in setdiff(setdiff(1:10,i),k)){mid = mid+nrow(filter(transmeta[!is.na(transmeta[,'d1']),],d1 == k & d5 ==h))}}
    for (k in setdiff(1:10,j)) {for(h in setdiff(setdiff(1:10,i),k)){mid = mid+nrow(filter(transmeta[!is.na(transmeta[,'d1']),],d1 == k & d5 ==h))}}
    
    trans[2,2] = mid
    #trans[2,2]<-nrow(filter(transmeta[!is.na(transmeta[,'d1']),],d1 != j & d5 != i))
    pva[10*(as.numeric(j)-1)+as.numeric(i),'trans_to_cs'] = as.numeric(i)
    pva[10*(as.numeric(j)-1)+as.numeric(i),'cs'] = as.numeric(j)
    pva[10*(as.numeric(j)-1)+as.numeric(i),'pvalue'] = fisher.test(trans)$p
    #print(paste('cluter',j,'trans to',i,'pvalue = ',fisher.test(trans)$p))
  }                
}

#pva$qvalue<-p.adjust(pva$pvalue,method = 'BH')

pva_sig<-filter(pva,pvalue <= 0.05)

pva_sig


# severity cs234 cs15789 --------------------------------------------------


##severity cs234 cs15789
metadata = read.csv('F:/ZLF/CAP/data/relative_data/metadata-all-220617.csv',row.names = 1)
trans<-data.frame(matrix(nrow = 2,ncol = 2))
trans<-data.frame(matrix(nrow = 2,ncol = 2))
rownames(trans)<-c('severe','non-severe')
colnames(trans)<-c('cs234','cs15789')

trans[1,1]<-filter(metadata,pam_10_cluster == 2 | pam_10_cluster ==3 | pam_10_cluster == 4) %>% 
  filter(severe_case == 'yes') %>% nrow()
trans[2,1]<-filter(metadata,pam_10_cluster == 2 | pam_10_cluster ==3 | pam_10_cluster == 4) %>% 
  filter(severe_case == 'no') %>% nrow()
trans[1,2]<-filter(metadata,pam_10_cluster == 1 | pam_10_cluster ==5 | pam_10_cluster == 7 | pam_10_cluster == 8 | pam_10_cluster == 9) %>% 
  filter(severe_case == 'yes') %>% nrow()
trans[2,2]<-filter(metadata,pam_10_cluster == 1 | pam_10_cluster ==5 | pam_10_cluster == 7 | pam_10_cluster == 8 | pam_10_cluster == 9) %>% 
  filter(severe_case == 'no') %>% nrow()
fisher.test(trans)


# severity cs234 cs6 ------------------------------------------------------


##severity cs234 cs6
metadata = read.csv('F:/ZLF/CAP/data/relative_data/metadata-all-220617.csv',row.names = 1)
trans<-data.frame(matrix(nrow = 2,ncol = 2))
trans<-data.frame(matrix(nrow = 2,ncol = 2))
rownames(trans)<-c('severe','non-severe')
colnames(trans)<-c('cs234','cs6')

trans[1,1]<-filter(metadata,pam_10_cluster == 2 | pam_10_cluster ==3 | pam_10_cluster == 4) %>% 
  filter(severe_case == 'yes') %>% nrow()
trans[2,1]<-filter(metadata,pam_10_cluster == 2 | pam_10_cluster ==3 | pam_10_cluster == 4) %>% 
  filter(severe_case == 'no') %>% nrow()
trans[1,2]<-filter(metadata,pam_10_cluster == 6) %>% 
  filter(severe_case == 'yes') %>% nrow()
trans[2,2]<-filter(metadata,pam_10_cluster == 6) %>% 
  filter(severe_case == 'no') %>% nrow()
fisher.test(trans)



# intubated patients were transmitted to CS5 in Day 5  --------------------


##intubated patients were transmitted to CS5 in Day 5 
metadata = read.csv('F:/ZLF/CAP/data/relative_data/metadata-all-220617.csv',row.names = 1)
trans<-data.frame(matrix(nrow = 2,ncol = 2))
trans<-data.frame(matrix(nrow = 2,ncol = 2))
rownames(trans)<-c('intubate','unintubate')
colnames(trans)<-c('tocs5','tocsothers')

transmeta_mv<-filter(transmeta,respiratory_support == '1' | respiratory_support == '2,1')
transmeta_oth<-filter(transmeta,respiratory_support != '1' & respiratory_support != '2,1')
trans[1,1]<-filter(transmeta,d5 == 5) %>% 
  filter(respiratory_support == '1' | respiratory_support == '2,1') %>% nrow()
trans[2,1]<-filter(transmeta,d5 == 5) %>% 
  filter(respiratory_support != '1' & respiratory_support != '2,1') %>% 
  nrow()
trans[1,2]<-filter(transmeta,d5 != 5) %>% 
  filter(respiratory_support == '1' | respiratory_support == '2,1') %>% nrow()
trans[2,2]<-filter(transmeta,d5 != 5) %>% 
  filter(respiratory_support != '1' & respiratory_support != '2,1') %>% 
  nrow()
fisher.test(trans)


