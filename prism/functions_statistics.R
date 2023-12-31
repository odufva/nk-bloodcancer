fisher.2x2=function(lv1, lv2, alternative="two.sided", usefast=T){
  cont=table(lv1, lv2)
  cont[is.na(cont)]=0
  
  if(usefast){
    p=HighSpeedStats::ultrafastfet(cont[1,1], cont[1,2], cont[2,1], cont[2,2])
    odds=NA
  }
  
  else{
    ft=fisher.test(cont, alternative=alternative)
    p=ft$p.value
    odds=ft$estimate
  }
  
  data.frame("p"=p, "log.odds"=odds, stringsAsFactors = F)
}

TestGeneWilcox=function(gene, datam, logicalVectors, logicalVector_normals=NULL, ALTERNATIVE="greater", prettynum=T){
  if(!gene%in%rownames(datam)) return(NULL)
  D = as.numeric(datam[rownames(datam)%in%gene,])
  
  d_list=lapply(logicalVectors, function(v){
    D[v][!is.na(D[v])]
  })
  
  # in this case, compute against rest of the samples for the logicalV
  if(is.null(logicalVector_normals)){
    logicalVector_normals=lapply(logicalVectors, '!') # negate
    names(logicalVector_normals)=rep("Rest", length(logicalVector_normals))
    
    d_list2=lapply(logicalVector_normals, function(v){
      D[v][!is.na(D[v])]
    })
    
    stats=do.call(rbind, lapply(seq(d_list), function(i, d_list, d_list2){
      l=d_list[[i]]
      name1=names(d_list)[i]
      name2=names(d_list2)[i]
      
      P.value=sapply(seq(d_list2)[i], function(j, l, d_list2){
        l2=d_list2[[j]]
        signif(wilcox.test(l, l2, ALTERNATIVE)$p.value,2)
      }, l, d_list2)
      
      FC=sapply(d_list2[i], function(l2, l1){
        2^(mean(l)-mean(l2))
      })
      
      df=data.frame("Gene"=gene, "Group1"=name1, "Group2"=name2, FC, "p"=P.value, "alternative"=ALTERNATIVE, stringsAsFactors=F)
      rownames(df)=NULL
      df=df[!df$Group1==df$Group2,]
      
      if(prettynum){
        df[,4]=prettyNum(signif(df[,4], 2))
        df[,5]=prettyNum(signif(df[,5], 2))
      }
      
      return(df)
    },d_list, d_list2))
    return(stats)
  }
  
  d_list2=lapply(logicalVector_normals, function(v){
    D[v][!is.na(D[v])]
  })
  
  stats=do.call(rbind, lapply(seq(d_list), function(i, d_list, d_list2){
    l=d_list[[i]]
    name1=names(d_list)[i]
    name2=names(d_list2)
    
    P.value=sapply(seq(d_list2), function(j, l, d_list2){
      l2=d_list2[[j]]
      signif(wilcox.test(l, l2, ALTERNATIVE)$p.value,2)
    }, l, d_list2)
    
    FC=sapply(d_list2, function(l2, l1){
      2^(mean(l)-mean(l2))
    })
    
    df=data.frame("Gene"=gene, "Group1"=name1, "Group2"=name2, FC, "p"=P.value, "alternative"=ALTERNATIVE, stringsAsFactors=F)
    rownames(df)=NULL
    df=df[!df$Group1==df$Group2,]
    
    if(prettynum){
      df[,4]=prettyNum(signif(df[,4], 2))
      df[,5]=prettyNum(signif(df[,5], 2))
    }
    
    return(df)
  },d_list, d_list2))
  
  return(stats)
}

FUN_HGT_MW_TEST=function(genes, logicalVector,logicalVector_normals, name, TEST_FAILS, ALTERNATIVE="greater", data = data, profile = profile, CORES=3, HG_PVAL=0.001, MW_PVAL=0.001, P.ADJUST.METH="bonferroni") {
  
  # transpose if needed
  if(any(grepl("TP53", rownames(data))))data=t(data)
  if(any(grepl("TP53", rownames(profile))))profile=t(profile)
  
  # Check some prerequisites
  if(nrow(data)!=nrow(profile)) stop("data and profile dimension mismatch")
  if(length(logicalVector)!=nrow(data)) stop("logicalVector dimension mismatch.")
  
  # Initialize stuff
  genes = genes[genes%in%colnames(data)]
  P = profile[,colnames(profile)%in%genes]
  D = data[,colnames(data)%in%genes]
  
  R = as.data.frame(matrix(NA, nrow = length(genes), ncol = 3))
  colnames(R) = c("gene", "pvalue", "adj.pvalue")
  R[,1] = genes
  
  
  #************ Hyperg test, disease association ***************
  FUN=function(j){
    r = genes[j]
    k = sum(logicalVector)
    x = P[,colnames(P)==r]
    m = sum(x)
    n = length(x) - m
    q = sum(x[logicalVector])
    phyper(q-1,m,n,k, lower.tail = F)
  }
  
  R[,2]=unlist(mclapply(1:nrow(R), FUN, mc.cores=CORES))
  
  R[,3] = p.adjust(R[,2], P.ADJUST.METH)
  
  # log10 transformation and filtering
  R = R[order(R$adj.pvalue),]
  R = R[R$adj.pvalue < HG_PVAL,] # FILTER
  R$pvalue = round(abs(log10(R$pvalue)),1)
  R$adj.pvalue = round(abs(log10(R$adj.pvalue)),1)
  
  #******************************************************************
  
  if(dim(R)[1]>0){
    
    # Mann-whitney, higher expression in normal
    MW = do.call(rbind, mclapply(1:dim(R)[1],function(i) {
      x=R[i,]
      gene = as.character(x[1])
      S = as.numeric(D[,colnames(D)==gene])
      logicalVector2=as.numeric(P[,colnames(D)==gene])
      
      # choose
      S_D = S[logicalVector&logicalVector2]
      S_D = S[logicalVector]
      
      ttp=sapply(logicalVector_normals, function(lv){
        S_normal = D[lv,colnames(D)==gene]
        tt = wilcox.test(S_D, S_normal,alternative=ALTERNATIVE)
        return(tt$p.value)
      })
      
      ttp = p.adjust(ttp, P.ADJUST.METH)
    }))
    
    testsig=MW>MW_PVAL
    
    name_normals=names(logicalVector_normals)
    
    res=do.call(cbind, lapply(seq(name_normals), function(i){
      ifelse(testsig[,i], name_normals[i], "")
    }))
    
    fails=apply(res, 1, function(v)paste(v[!v==""], collapse=","))
    
    # also compute fold change against all normals
    options("scipen"=100, "digits"=3)
    FC = do.call(rbind, mclapply(1:dim(R)[1],function(i) {
      x=R[i,]
      gene = as.character(x[1])
      S = as.numeric(D[,colnames(D)==gene])
      logicalVector2=as.numeric(P[,colnames(D)==gene])
      
      # choose
      S_D = S[logicalVector&logicalVector2]
      S_D = S[logicalVector]
      
      fc=sapply(logicalVector_normals, function(lv){
        S_normal = D[lv,colnames(D)==gene]
        FC=mean(S_normal)-mean(S_D)
      })
    }, mc.cores=CORES))
    
    FoldChange=round(rowMeans(FC)*-1,1)
    
    FoldChange_min=round(apply(FC, 1, min)*-1,1)
    FoldChange_max=round(apply(FC, 1, max)*-1,1)
    
    # report all of these
    df=data.frame(R, FoldChange, fails, FoldChange_min, FoldChange_max)
    
  }else{
    return(NULL)
  }
  
  # output
  if(dim(df)[1]>0)  return(cbind(df, name, stringsAsFactors=F))
}

fun.kapplanMeier=function(TIME, STATUS, GROUPS=NULL, CONTINUOUS=NULL, CONTINUOUS_SUMMARY=c("maxrank", "5quantiles","4quantiles","3quantiles","85th_15th_percentile", "80th_20th_percentile", "75th_25th_percentile", "median"), NAME="", MONTHS=F, PVAL=1, INDIVIDUAL_GROUPS=T, LWD=3, labels=NULL, COLORV=NULL){
  
  CONTINUOUS_SUMMARY <- match.arg(CONTINUOUS_SUMMARY)
  
  # cutpoints https://github.com/kassambara/survminer/issues/41
  if(MONTHS){
    TIME=as.numeric(TIME)*0.0328767
  }
  
  if(is.null(GROUPS)&!is.null(CONTINUOUS)){
    df=data.frame("time"=TIME, "status"=STATUS, "continuous"=CONTINUOUS, stringsAsFactors = F)
    
    NAME=paste(NAME, CONTINUOUS_SUMMARY)
    
    if(CONTINUOUS_SUMMARY=="maxrank"){
      library(survminer)
      cutoff <- surv_cutpoint(
        df,
        time = "time",
        event = "status",
        variables = c("continuous")
      )
      # print(plot(cutoff))
      GROUPS=rep(paste0("low (<", signif(cutoff$cutpoint$cutpoint,1), ")"), length(CONTINUOUS))
      GROUPS[CONTINUOUS>cutoff$cutpoint$cutpoint]=paste0("High (>=", signif(cutoff$cutpoint$cutpoint,1), ")")
    }
    
    if(CONTINUOUS_SUMMARY=="4quantiles"){
      
      library(dplyr)
      GROUPS <- paste0("Q",ntile(CONTINUOUS, 4))
      
    }
    
    if(CONTINUOUS_SUMMARY=="3quantiles"){
      
      library(dplyr)
      GROUPS <- paste0("Q",ntile(CONTINUOUS, 3))
      
    }
    
    if(CONTINUOUS_SUMMARY=="5quantiles"){
      
      library(dplyr)
      GROUPS <- paste0("Q",ntile(CONTINUOUS, 5))
      
    }
    
    
    if(CONTINUOUS_SUMMARY=="85th_15th_percentile"){
      
      q=quantile(CONTINUOUS, probs = c(0.15, 0.85)) # quartile
      GROUPS=rep("Q2", length(CONTINUOUS))
      GROUPS[CONTINUOUS>q[2]]="Q3"
      GROUPS[CONTINUOUS<q[1]]="Q1"
      
    }
    
    if(CONTINUOUS_SUMMARY=="80th_20th_percentile"){
      
      q=quantile(CONTINUOUS, probs = c(0.2, 0.8)) # quartile
      GROUPS=rep("Q2", length(CONTINUOUS))
      GROUPS[CONTINUOUS>q[2]]="Q3"
      GROUPS[CONTINUOUS<q[1]]="Q1"
      
    }
    
    if(CONTINUOUS_SUMMARY=="75th_25th_percentile"){
      
      q=quantile(CONTINUOUS, probs = c(0.25, 0.75)) # quartile
      GROUPS=rep("Q2", length(CONTINUOUS))
      GROUPS[CONTINUOUS>q[2]]="Q3"
      GROUPS[CONTINUOUS<q[1]]="Q1"
      
    }
    
    if(CONTINUOUS_SUMMARY=="median"){
      
      library(dplyr)
      GROUPS <- paste0("Q",ntile(CONTINUOUS, 2))
      
    }
  }
  
  classes=as.character(GROUPS)
  size=c(length(unique(classes)))
  
  if(!is.null(COLORV)){
    col1 <- COLORV
  }else{
    size=c(length(unique(classes)))
    col1=c("red", "gray75")
    if(size>2)col1 <- colorRampPalette(c(brewer.pal(min(size, 9), "Set1"), brewer.pal(min(size, 8), "Set2"), brewer.pal(min(size, 9), "Set3")))(size)
  }
  
  if(!INDIVIDUAL_GROUPS){
    DF.surv2=data.frame("time"=as.numeric(TIME), "status"=as.character(STATUS)=="DECEASED"|as.character(STATUS)=="1"|as.character(STATUS)=="Dead", "x"=as.character(GROUPS), stringsAsFactors=F)
    d.surv2 <- survfit(Surv(time,status) ~ x, data = DF.surv2, type="kaplan-meier", conf.type="log")
    
    plot(d.surv2, col=col1, lwd=rep(as.numeric(LWD), length(col1)), mark.time=c(as.numeric(TIME)), mark=3, cex=LWD/2, fun="log", conf.int=F, ylim=c(0, 1), ylab="Percentage Surviving", xlab="Months")
    legend("topright", paste(names(d.surv2$strata)," n=",d.surv2$n, sep=""), lwd=rep(LWD, length(col1)), col=col1)
    title(paste("Kaplan-Meier Curves\n", NAME))
  }else{
    
    unfeat=unique(as.character(GROUPS))
    # unfeat=unfeat[!unfeat==F&!grepl("0_", unfeat)]
    
    for(a in unfeat){
      print(a)
      group=ifelse(GROUPS%in%a, a, "Rest")
      
      
      DF.surv2=data.frame("time"=as.numeric(TIME), "status"=as.character(STATUS)=="DECEASED"|as.character(STATUS)=="1"|as.character(STATUS)=="Dead", "x"=as.character(group), stringsAsFactors=F)
      d.surv2 <- survfit(Surv(time,status) ~ x, data = DF.surv2, type="kaplan-meier", conf.type="log")
      
      classes=as.character(GROUPS)
      size=c(length(unique(classes)))
      
      if(!is.null(COLORV)){
        col1 <- COLORV
      }else{
        size=c(length(unique(classes)))
        col1=c("red", "gray75")
      }
      
      # why is this here?
      # if(all(table(DF.surv2[,2], DF.surv2[,3])>=2)){
      Hazard.t=coxph(Surv(time,status) ~ x, data = DF.surv2)
      add=signif(summary(Hazard.t)$logtest["pvalue"],2)[1]
      # }else{
      #   add=NA
      # }
      
      if(!is.na(add)&add<=PVAL){
        plot(d.surv2, col=col1, lwd=rep(as.numeric(LWD), length(col1)), mark.time=c(as.numeric(TIME)), mark=3,cex=LWD/2, fun="log", conf.int=F, ylim=c(0, 1), xlab="Months", ylab="Percentage Surviving")
        legend("topright", paste(gsub("x=", "", names(d.surv2$strata))," n=",d.surv2$n, " ", c("", paste0("p.value=", add)), sep=""), lwd=rep(LWD, length(col1)), col=col1)
        title(paste("Kaplan-Meier Curves\n", NAME))
        
      }
    }
  }
}
fun.coxph=function(TIME, STATUS, GROUPS, NAME="", MONTHS=F, PVAL=1, INDIVIDUAL_GROUPS=T, labels=NULL, COLORV=NULL){
  if(MONTHS){
    TIME=as.numeric(TIME)*0.0328767
  }
  
  library("survival")
  library("survminer")
  
  res.cox <- coxph(Surv(TIME, STATUS) ~ age + sex + ph.ecog, data =  lung)
  summary(res.cox)
  
  # Survival curves
  fit <- survfit(res.cox, newdata = sex_df)
  ggsurvplot(fit, conf.int = TRUE, legend.labs=c("Sex=1", "Sex=2"),
             ggtheme = theme_minimal())
  
  classes=as.character(GROUPS)
  size=c(length(unique(classes)))
  
  if(!is.null(COLORV)){
    col1 <- COLORV
  }else{
    size=c(length(unique(classes)))
    col1=c("red", "gray75")
    if(size>2)col1 <- colorRampPalette(c(brewer.pal(min(size, 9), "Set1"), brewer.pal(min(size, 8), "Set2"), brewer.pal(min(size, 9), "Set3")))(size)
  }
  
  if(!INDIVIDUAL_GROUPS){
    DF.surv2=data.frame("time"=as.numeric(TIME), "status"=as.character(STATUS)=="DECEASED"|as.character(STATUS)=="1"|as.character(STATUS)=="Dead", "x"=as.character(GROUPS), stringsAsFactors=F)
    d.surv2 <- survfit(Surv(time,status) ~ x, data = DF.surv2, type="kaplan-meier", conf.type="log")
    
    plot(d.surv2, col=col1, lwd=rep(3, length(col1)), mark.time=c(as.numeric(TIME)), mark=3, fun="log", conf.int=F, ylim=c(0, 1), ylab="Percentage Surviving", xlab="Months")
    legend("topright", paste(names(d.surv2$strata)," n=",d.surv2$n, sep=""), lwd=rep(4, length(col1)), col=col1)
    title(paste("Kaplan-Meier Curves", NAME))
  }else{
    
    unfeat=unique(as.character(GROUPS))
    unfeat=unfeat[!unfeat==F&!grepl("0_", unfeat)]
    
    for(a in unfeat){
      print(a)
      group=ifelse(GROUPS%in%a, a, "Rest")
      
      
      DF.surv2=data.frame("time"=as.numeric(TIME), "status"=as.character(STATUS)=="DECEASED"|as.character(STATUS)=="1"|as.character(STATUS)=="Dead", "x"=as.character(group), stringsAsFactors=F)
      d.surv2 <- survfit(Surv(time,status) ~ x, data = DF.surv2, type="kaplan-meier", conf.type="log")
      
      classes=as.character(GROUPS)
      size=c(length(unique(classes)))
      
      if(!is.null(COLORV)){
        col1 <- COLORV
      }else{
        size=c(length(unique(classes)))
        col1=c("red", "gray75")
      }
      
      if(all(table(DF.surv2[,2], DF.surv2[,3])>=2)){
        Hazard.t=coxph(Surv(time,status) ~ x, data = DF.surv2)
        add=signif(summary(Hazard.t)$logtest["pvalue"],2)[1]
      }else{
        add=NA
      }
      
      if(!is.na(add)&add<=PVAL){
        plot(d.surv2, col=col1, lwd=rep(3, length(col1)), mark.time=c(as.numeric(TIME)), mark=3, fun="log", conf.int=F, ylim=c(0, 1), xlab="Months", ylab="Percentage Surviving")
        legend("topright", paste(gsub("x=", "", names(d.surv2$strata))," n=",d.surv2$n, " ", c("", paste0("p.value=", add)), sep=""), lwd=rep(4, length(col1)), col=col1)
        title(paste("Kaplan-Meier Curves", NAME))
        
      }
    }
  }
}

cor.test.p=function(x, y, method="spearman", use="pairwise.complete.obs") cor.test(x, y, method=method, use=use)[["p.value"]]

cor.test.p.m <- function(x, method="spearman", use="pairwise.complete.obs"){
  
  z <- outer(
    colnames(x),
    colnames(x),
    Vectorize(function(i,j){
      no.pairs=sum(complete.cases(x[,j]==x[,i]))<3
      if(no.pairs)return(NA)
      cor.test.p(x[,i], x[,j], method=method, use=use)
    })
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

flat_cor_mat <- function(cor_rho, cor_pval){
  # cor_3 <- rcorr(as.matrix(mtcars[, 1:7]))
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(dplyr)
  library(tibble)
  options(scipen=4)
  cor_rho <- rownames_to_column(as.data.frame(cor_rho), var = "row")
  cor_rho <- gather(cor_rho, column, cor, -1)
  cor_pval <- rownames_to_column(as.data.frame(cor_pval), var = "row")
  cor_pval <- gather(cor_pval, column, p, -1)
  cor_pval_matrix <- left_join(cor_rho, cor_pval, by = c("row", "column"))
  colnames(cor_pval_matrix)[1:2]=c("featureA", "featureB")
  cor_pval_matrix=cor_pval_matrix[!cor_pval_matrix[,1]==cor_pval_matrix[,2],]
  cor_pval_matrix=cor_pval_matrix[order(cor_pval_matrix[,1], cor_pval_matrix[,3]),]
  rownames(cor_pval_matrix)=NULL
  cor_pval_matrix
}

chart.Correlation.custom=function (R,filterv=NULL, histogram = TRUE, method = c("pearson", "kendall","spearman"), ...){
  x = checkData(R, method = "matrix")
  if (missing(method))
    method = method[1]
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs",
                        method, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor))
      cex <- 0.8/strwidth(txt)
    test <- cor.test(x, y, method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***",
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  hist.panel = function(x, ...) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, axes = FALSE,
         main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  if (histogram)
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,
          diag.panel = hist.panel, method = method, subset=colnames(x)%in%rownames(test)[1:5])
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,
             method = method, ...)
}
