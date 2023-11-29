regulon.feats <- function(fm, genelist, cnv_annot=NULL, filter.val=3, filtertypes=NULL){
  fmname=rownames(fm)
  
  if(!is.null(filtertypes)){
    fmname=fmname[!substr(fmname, 1,6)%in%filtertypes]
  }
  feature.gene=suppressWarnings(do.call(rbind, strsplit(fmname, ":"))[,3])
  
  if(any(grepl(":", genelist)))genelist[grepl(":", genelist)]=do.call(rbind, strsplit(genelist[grepl(":", genelist)], ":"))[,3]
  
  if(!is.null(cnv_annot)){
    l.genes=strsplit(cnv_annot[,2], ",")
    names(l.genes)=cnv_annot[,1]
    st=stack(l.genes)
    st=st[st[,1]%in%genelist,]
    st[,2]=as.character(st[,2])
  }
  
  nested.features=function(gene, fmname, feature.gene){
    direct=fmname[feature.gene%in%gene]
    if(!is.null(cnv_annot))direct=c(direct, st[st[,1]%in%gene,2])
    
    # order this to numeric, then to gexp-cnv-meth-gnab
    n=unique(direct)
    types=substr(n, 1, 6)
    ord=c("N:GEXP","B:GEXP","B:GNAB", "N:CNVR", "B:CNVR","N:METH", "N:MIRN", "N:SAMP", "B:SAMP", "N:CLIN", "B:CLIN", "N:DRUG", "N:PRSM")
    ord=ord[ord%in%types]
    
    r=unlist(lapply(ord, function(o) n[types%in%o]))
    
    return(r)
  }
  feats=lapply(genelist, nested.features, fmname, feature.gene)
  
  # Filter if not enough data:
  m=fm[rownames(fm)%in%unlist(feats),]
  s=is.na(m)
  m=m[!rowSums(s)==dim(m)[1]&rowSums(!s)>=filter.val, !colSums(s)==dim(m)[2]]
  
  feats=lapply(feats, function(n)n[n%in%rownames(m)])
  names(feats)=genelist
  
  feats=feats[!unlist(lapply(feats, length))<1]
  
  return(feats)
}
