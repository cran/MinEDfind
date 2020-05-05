select.MinED <- function(n, y, z, phi_t, phi_e, eps_t, eps_e, ct = 0.95){
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
  # only dose level which has the trial will be considered into the admissable set
  n_trial <- n[n > 0]
  admissable_set <- admtox(n_trial, y, ct, phi_t)
  dselect <- NA
  if (admissable_set[1] == 0){
    print("No dose selected")
  }
  else{
    c_maxid <- rep(NA, length(n))
    target_e <- rep(NA, length(n))
    post_efficacy_dist_generator <- vector("list", length(n))
    post_toxicity_dist_generator <- vector("list", length(n))
    posterior_efficacy_mean <- rep(NA, length(n))
    posterior_toxicity_mean <- rep(NA, length(n))
    upper_eff = rep(NA, length(n))
    lower_eff = rep(NA, length(n))
    upper_tox = rep(NA, length(n))
    lower_tox = rep(NA, length(n))
    for(id in 1:length(n)){
      post_efficacy_dist_generator[[id]] <- function(q){
        q * dbeta(q, z[id]+1, n[id]-z[id]+1)
      }
      post_toxicity_dist_generator[[id]] <- function(w){
        w * dbeta(w, y[id]+1, n[id]-y[id]+1)
      }
      c_maxid[id] = maxKey(n[id],z[id],phi_e,eps = eps_e)$idmax
      target_e[id] = maxKey(n[id],z[id],phi_e,eps = eps_e)$target
      posterior_efficacy_mean[id] = integrate(post_efficacy_dist_generator[[id]], 0, 1)$value
      posterior_toxicity_mean[id] = integrate(post_toxicity_dist_generator[[id]], 0, 1)$value
      upper_eff[id] = qbeta(0.975, z[id]+1, n[id]-z[id]+1)
      lower_eff[id] = qbeta(0.025, z[id]+1, n[id]-z[id]+1)
      upper_tox[id] = qbeta(0.975, y[id]+1, n[id]-y[id]+1)
      lower_tox[id] = qbeta(0.025, y[id]+1, n[id]-y[id]+1)
    }
    dis = (c_maxid[1:admissable_set[2]]-target_e[1:admissable_set[2]])
    if(any(dis>=0)){
      dselect = which(dis>=0)[1]
    }
    else{
      dselect= which.min(abs(dis))[1]
    }
  }
  # dose level which doesn't have trials won't have the information below
  posterior_efficacy_mean[n == 0] <- NA
  lower_eff[n == 0] <- NA
  upper_eff[n == 0] <- NA
  posterior_toxicity_mean[n == 0] <- NA
  lower_tox[n == 0] <- NA
  upper_tox[n == 0] <- NA
  target_e[n == 0] <- NA
  c_maxid[n == 0] <- NA
  df_info <- data.frame(Dose_Level = 1:length(n), Posterior_Efficacy_Est =  unlist(posterior_efficacy_mean),
                        Lower_Efficacy = lower_eff, Upper_Efficacy = upper_eff, Posterior_Toxicity_Est = unlist(posterior_toxicity_mean),
                        Lower_Toxicity = lower_tox, Upper_Toxicity = upper_tox, Target_ID = target_e, MaxKey_ID = c_maxid)
  target_level <- data.frame(Target_Efficacy_Rate = phi_e, Target_Toxicity_Rate = phi_t)
  outcome <- list(Selected_Dose = dselect, Target_Level = target_level, Info = df_info)
  class(outcome) <- c("MinED", "list")
  outcome
}
