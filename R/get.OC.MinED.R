get.OC.MinED <- function(ttox, teff, phi_t, phi_e, ct=0.95, eps_t, eps_e, d0 = 1,
                        cohortsize = 3, ncohort1, ncohort2, ntrial = 100, extrasafe = TRUE,
                        cutoff.eli = 0.95, n.earlystop = 12){
  maxKey <- function(n,y,phi,eps){
    prkeyboard = NULL
    cutpoint = sort(c(seq(phi-eps,0,by=-eps * 2),seq(phi+eps,1,by=eps * 2)),decreasing = F) ## by = eps * 2
    len = length(cutpoint)
    key_L = cutpoint[-len]
    key_R = cutpoint[-1]
    for(i in 1:(len-1)){
      prkeyboard[i] = pbeta(key_R[i],y+1,n-y+1) - pbeta(key_L[i],y+1,n-y+1)
    }
    #idmax= which.max(prkeyboard)
    idmax= which.max(round(prkeyboard, 10)) ## make sure that R will ignore very small difference due to the computation limit
    id_target = which(key_L==phi-eps)
    list(idmax=idmax, target=id_target, v_max=prkeyboard[idmax])
  }
  admtox <- function(n,y,ct,phi_t){
    dtried = max(which(n!=0))
    rest = rep(0, length(n))
    phi_t = phi_t#+0.1*phi_t
    rest[1:dtried]=1-pbeta(phi_t,1+y[1:dtried],1+n[1:dtried]-y[1:dtried])
    rest=pava(rest,w=n)
    if(any(rest<=ct)){
      drange=range(which(rest<=ct))
    }else{
      drange = 0
    }
    return(drange)
  }

  keylead <- function(d,drange,peak.d,dtried,max_e,target_e){
    if(max_e<target_e){
      if(d<peak.d){
        d= min(d+1,drange[2])
      }else{
        if(d>peak.d){
          d= max(d-1,drange[1])
        }else{
          if(peak.d==dtried){
            d=min(d+1,drange[2])
          }else{
            d=d
          }
        }
      }
    }else{
      if(max_e>target_e){
        d = max(d-1, drange[1])
      }else{
        d = d
      }
    }
    return(d)
  }

  ncohort = ncohort1+ncohort2
  ndose=length(ttox)
  N=matrix(rep(0,ndose*ntrial),ncol=ndose)
  YTOX=matrix(rep(0,ndose*ntrial),ncol=ndose)
  YEFF=matrix(rep(0,ndose*ntrial),ncol=ndose)
  dselect = rep(0, ntrial)
  for(itrial in 1:ntrial){
    flag = "green"
    y=rep(0,ndose)
    z=rep(0,ndose)
    n=rep(0,ndose)
    d = d0
    for (i in 1:ncohort){
      if(i<=ncohort1){
        y[d]=y[d]+rbinom(1,cohortsize,ttox[d])
        z[d]=z[d]+rbinom(1,cohortsize,teff[d])
        n[d]=n[d]+cohortsize
        if(1-pbeta(phi_t,1+y[1],1+n[1]-y[1])>cutoff.eli&n[1]>n.earlystop&extrasafe){
          dselect[itrial] = 99
          flag="red"
          break
        }
        resKey_t = maxKey(n[d],y[d],phi=phi_t,eps=eps_t)
        max = resKey_t$idmax
        target = resKey_t$target
        if(max<target){d = min(d+1, ndose)}
        if(max>target) d = max(d-1,1)
        if(i==ncohort1){
          drange= admtox(n,y,ct,phi_t)
          if(drange[1]==0){
            dselect[itrial] = 99
            flag="red"
            break
          }else{
            c_d = drange[1]:drange[2]
            c_maxid = NULL
            for(id in 1:length(c_d)){
              c_maxid[id] = maxKey(n[id],z[id],phi_e,eps = eps_e)$idmax
              target_e = maxKey(n[id],z[id],phi_e,eps = eps_e)$target
            }
            dis = (c_maxid-target_e)
            if(any(dis>=0)){
              d = which(dis>=0)[1]
            }else{
              d = which.min(abs(dis))[1]
            }
          }
          }
      }else{
        y[d]=y[d]+rbinom(1,cohortsize,ttox[d])
        z[d]=z[d]+rbinom(1,cohortsize,teff[d])
        n[d]=n[d]+cohortsize
        drange= admtox(n,y,ct,phi_t)
        if(1-pbeta(phi_t,1+y[1],1+n[1]-y[1])>0.95&n[1]>6&extrasafe){
          dselect[itrial] = 99
          flag="red"
          break
        }
        resKey_e = maxKey(n[d],z[d],phi=phi_e,eps=eps_e)
        max_e = resKey_e$idmax
        target_e = resKey_e$target
        if(drange[1]!=0){
          dtried=max(which(n>0))
          zfit = (z+0.01)/(n+0.01)
          unifit=ufit(zfit, x=1:ndose, w=n, type="b")
          dopt = unifit$mode
          d =  keylead(d,drange,dopt,dtried,max_e,target_e)
        }else{
          dselect[itrial] = 99
          flag="red"
          break
        }
      }
    }
    if(flag=="green"){
          c_maxid = NULL
          dtried = max(which(n>0))
          for(id in 1:dtried){
            c_maxid[id] = maxKey(n[id],z[id],phi_e,eps = eps_e)$idmax
          }
          dis = (c_maxid-target_e)
          if(any(dis>=0)){
            dselect[itrial] = which(dis>=0)[1]
            }else{
              dselect[itrial] = which.min(abs(dis))[1]
            }
    }
    N[itrial,]=n
    YTOX[itrial,]=y
    YEFF[itrial,]=z
  }
  selpercent=rep(0, ndose)
  for(i in 1:ndose) { selpercent[i]=sum(dselect==i)/ntrial*100 }
  oc<- matrix(NA,ncol = ndose+1, nrow = 6)
  oc[1,1:ndose] = ttox
  oc[2,1:ndose] = teff
  oc[3,1:ndose] = selpercent
  oc[4,1:(ndose+1)] = c(apply(N,2,mean),sum(apply(N,2,mean)))
  oc[5,1:(ndose+1)] = c(apply(YTOX,2,mean),sum(apply(YTOX,2,mean)))
  oc[6,1:(ndose+1)] = c(apply(YEFF,2,mean),sum(apply(YEFF,2,mean)))
  colnames(oc)=c(paste("Dose",1:ndose,sep=""),"Total # Pts")
  row.names(oc) =c("True DLT rate","True efficacy rate","Sec %",
                   "# Pts treated","# Pts response to tox","# Pts response to eff")
  class(oc) <- c("MinED", "matrix")
  return(oc)
}
