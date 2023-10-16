# Function to calculate Shannon-entropy of vector
shannonEntropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm<-p[p>0]/sum(p)
  -sum( log2(p.norm)*p.norm)
}

# Function to convert FPKM/counts into proportions
makeProbs<-function(a){
  colSums<-apply(a,2,sum)
  b<-t(t(a)/colSums)
  b[is.na(b)] = 0
  b
}

# Function to calculate Jensen-Shannon distance
JSdist<-function(mat,...){
  res<-matrix(0,ncol=dim(mat)[2],nrow=dim(mat)[2])
  
  col_js<-apply(mat,MARGIN=2,shannonEntropy)
  
  colnames(res)<-colnames(mat)
  rownames(res)<-colnames(mat)
  for(i in 1:dim(mat)[2]){
    for(j in i:dim(mat)[2]){
      a<-mat[,i]
      b<-mat[,j]
      JSdiv<-shannonEntropy((a+b)/2)-(col_js[i]+col_js[j])*0.5
      res[i,j] = sqrt(JSdiv)
      res[j,i] = sqrt(JSdiv)
    }
  }
  res<-as.dist(res,...)
  attr(res,"method")<-"JSdist"
  res
}
# Function that does the clustering
# m - a matrix of FPKM (or TPM, or whatever)
# k - number of desired clusters
# logMode - whether or not to use the log10 of FPKM instead of raw FPKM values
# method - an alternative function to calculate distance between genes (default is Jensen-Shannon distance)
# pseudocount - this value will be added to every FPKM (to account for cases where FPKM = 0)
clusterExpr <- function(m, k, logMode=T, method='none', pseudocount=1,...){
  if(!require(cluster)) stop('Please install the "cluster" package with install.package("cluster")')
  
  m <- m[rowSums(m) > 0,]
  if(logMode){
    m<-log10(m+pseudocount)
  }
  
  if(!is.function(method)){
    method = function(mat){JSdist(makeProbs(t(m)))}	
  }		
  
  n<-method(m)
  
  clusters<-pam(n,k, ...)
  
  class(clusters)<-"list"
  clusters$fpkm<-m
  clusters
}
