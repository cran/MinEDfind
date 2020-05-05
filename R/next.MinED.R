next.MinED<-function(n, y, z, d, phi_t=0.3, phi_e = 0.3, eps_t, eps_e, ct=0.95, N1 =18){
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
      }else{# d >= peak.d
        if(d>peak.d){
          d= max(d-1,drange[1])
        }else{#d=peak.d
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

  ndose=length(y)
  adm=admtox(n,y,ct,phi_t)
  if(sum(n)<N1){
    resKey_t = maxKey(n[d],y[d],phi=phi_t,eps=eps_t)
    max = resKey_t$idmax
    target = resKey_t$target
    if(max<target){d = min(d+1, ndose)}
    if(max>target) d = max(d-1,1)
  }else{
    if(sum(n)==N1){
      drange= admtox(n,y,ct,phi_t)
      if(drange[1]==0){
        d <- 99
        flag="red"
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
    }else{
      drange= admtox(n,y,ct,phi_t)
      resKey_e = maxKey(n[d],z[d],phi=phi_e,eps=eps_e)
      max_e = resKey_e$idmax
      target_e = resKey_e$target
      if(drange[1]!=0){
        dtried=max(which(n>0))
        zfit = (z+0.01)/(n+0.01)
        unifit=ufit(zfit, x=1:ndose, w=n, type="b")
        peak.d = unifit$mode
        d = keylead(d, drange, peak.d, dtried, max_e, target_e)
      }else{
        d = 99
      }
    }
  }
  return(list("nextdose"=d))
}
