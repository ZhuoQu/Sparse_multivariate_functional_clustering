
generate_cluster_minprop<-function(nbvar, n, nbtime, clusternumber, 
                                   variat, partitionidentical) {
  part <- rep(NA, nbvar)
  if (partitionidentical == 1) {
    part[1] <- round(n * runif(1, 0.05, 0.8))
    part[2] <- round(n * runif(1, 0.05, 1-part[1]/n))
    part[3] <- n - part[1] - part[2]
  } else if (partitionidentical == 2){
    part[1] <- round(n * runif(1, 0.05, 0.45))
    part[2] <- part[1]
    part[3] <- n-2 * part[1]
  } else if (partitionidentical == 3) {
    part[1] <- round(n/clusternumber)
    part[2] <- part[1]
    part[3] <- n - part[1] - part[2]
  }
  
  var<-list()
  ### set  clusters
  for (k in 1:nbvar) {
    l <- (sample(2: nbvar,1) + sample(1: nbvar, 1))/2
    ran <- ifelse(rbinom(1, 1, 0.5) < 1, -1, 1)
    mat <- c()
    
    for (i in 1:clusternumber) { 
      
      if (variat == "shift") {
        
        val <- 6* cos((l + k) / 2 * pi * (t+0.2*ran*i))
      } else if (variat == "amplitude")
      {
        r<-ifelse(rbinom(1, 1, 0.5) < 1, -1, 1)
        val <- 2 * k * cos((l + ran * k/4) * pi * t)+(-1)^r * 3 * i * k
        
      }else if (variat == "phase")
      {
        if (i==1)
        {
          ht = log(t+1)+1-log(2)
        } else if (k==1 && i==2)
        {
          ht = t^2
        } else if (k==2 && i==2)
          {
          ht = 1 - cos(pi*t/2)
          } else if (k == 3 && i == 2)
            {
            ht = (sin(pi*t/2))^2
            }else if ( k == 1 && i == 3 ) {
          ht = t^3
        } else if ( k == 2 && i == 3 )
        {
          ht = sin (pi*t/2)
        } else if ( k == 3 && i == 3 )
        {
          ht = t
        }
        val<-2 * k * cos((l + k/4) * pi * ht)
        
      }else if (variat == "clover")
      {
        Theta<-seq(0,(2 * pi),length.out = 150)
        theta<-Theta[(14: 38) + (i-1) * 25]
        if (k == 1)
        {
          val<-5 * cos(3 * theta) * cos(theta)
        } else if (k == 2){
          val<-5 * cos(3 * theta) * sin(theta)
        } else
        {
          val<-5 * cos(3 * theta)
        }
      } else if (variat == "cyclone")
      {
        if (k==1 && i==1)
        {
          val<- -5 * t-sin(t * (10 * pi))
        } else if (k == 1 && i == 2)
          {
          val<-5 * sin(t * (10 * pi))
          } else if (k==1 &i==3)
          {
            val<-5*t+sin(t*(5*2*pi))
          }
        else if ( k==2 && i==1){
          val<-11 * t^2
        } else if (k==2 && i==2)
        {
          val<-7 * t
        } else if (k==2 && i==3)
        {
          val<-5 * log(t+1)
        }
        else if (k == 3) {
          val<- 7 * i * log(t+1)
        }
      } else if (variat == "helix") {
        theta <- 10+t*(10*2*pi)
        r <- 5*t
        if (k==1 && i==1)
        {
          val <- 5.5*cos(t*(5*2*pi))
        } else if ((k == 1 && i == 2)||(k == 1 && i == 3))
        {
          val <- r*cos(theta)
        } else if (k == 2 && i == 1)
        {
          val <- 5.5*sin(t*(5*2*pi))
        } else if ( (k==2 && i==2) ||(k==2 & i==3)){
          val <- r*sin(theta)
          
        } else if (k==3 &i==2)
        {
          val<- 10-2*r
        }
        else if (k==3){
          val <- 2*r
        }
      }
      
      matrix_element <- c()
      for (nn in 1: part[i]) {
        matrix_element <- c(matrix_element, runif(1, 0.9, 1.1) * val)
      }
      sub<-matrix(matrix_element, ncol = nbtime, nrow = part[i], byrow = TRUE)
      mat<-rbind(mat,sub)
    }
    var[[k]]<-mat
  }
  
  ########################################### 
  nu<-nu_construct(0.2,0.3,nbvar)
  beta<-beta_construct(nbvar)
  rho <- rho_construct(beta,nu,nbvar)
  a<-a_construct(2,nbvar)
  
  
  Sigma<-cross_covariance(nbvar,t,var=runif(nbvar,0.1,0.2),rho,nu,a)
  
  for (i in 1:n)
  {
    error<-mvrnorm(1,rep(0,nbvar*nbtime),Sigma)
    for (j in 1:nbvar)
    {
      var[[j]][i,]<-var[[j]][i,]+error[((j-1)*nbtime+1):(j*nbtime)]
    }
  }
  label <- c(rep(1, part[1]),rep(2, part[2]),rep(3, part[3]))
  return (list(value=var,label))
}

introduce_outlier<-function(var, outlier_nb, contamination_type){
  summary_min<-lapply(var[[1]], function(k){min(k[k!=min(k)])})
  summary_max<-lapply(var[[1]],function(k){max(k[k!=max(k)])})
  ##### line abnormality
  n <- nrow(var[[1]][[1]])
  
  if (length(outlier_nb)>0)
  {
    for (k in 1:nbvar)
    {
      rgd<-summary_max[[k]]-summary_min[[k]]
      if (contamination_type=="shape outlier I")
      {
        for (curve in outlier_nb)
        {
          l<-runif(1,-min(abs(summary_max[[k]]),abs(summary_min[[k]])),min(abs(summary_max[[k]]),abs(summary_min[[k]])))
          if (k==1) {
          var[[1]][[k]][curve,]<-var[[1]][[k]][curve,]+2*cos(2*pi*t) ## shape outlier I
          } else if (k==2) {
            var[[1]][[k]][curve,]<-var[[1]][[k]][curve,]+2*sin(2*pi*t) ## shape outlier I
          } else if (k==3) {
            var[[1]][[k]][curve,]<-var[[1]][[k]][curve,]+2*t ## shape outlier I
            
          }
          }
      } else if (contamination_type=="shape outlier II")
      {
        
        for (curve in outlier_nb)
        {
          constant<-runif(1,summary_min[[k]],summary_max[[k]]/2)
          slope<-runif(1,-2*k,2*k)
          var[[1]][[k]][curve,]<-constant+slope*t+runif(nbtime,-0.3,0.3) ## shape outlier II
        }
      } else if (contamination_type=="shape outlier III")
      {
        for (curve in outlier_nb)
        {
          l<-runif(1,-0.5,0.5)
          var[[1]][[k]][curve,]<-var[[1]][[k]][curve,]+((1+l)*k/2)*cos(10*nbvar*pi*t)
        }
      } else if (contamination_type=="pure magnitude outlier")
      {
        for (curve in outlier_nb)
        {
          shift_val<-ifelse(runif(1,-1,1)>0,1,-1)
          var[[1]][[k]][curve,]<-var[[1]][[k]][curve,]+shift_val*rep(rgd,nbtime)
        }
      }
      #  else if (contamination_type=="weakly pure magnitude outlier")
      # {
      #   shift_val<-ifelse(runif(1,-1,1)>0,1/3,-1/3)
      #   for (curve in outlier_nb)
      #   {
      #     var[[k]][curve,]<-var[[k]][curve,]+shift_val*rep(rgd,nbtime)
      #   }
      # }
      else if (contamination_type=="peak magnitude outlier")
      {
        for (curve in outlier_nb)
        {
          shift_val<-ifelse(runif(1,-1,1)>0,1/2,-1/2)
          rv<-sample(1:length(t),1)
          st<-ifelse(rv>0.9*length(t),floor(0.9*length(t))-2,rv)
          tp<-st:(st+0.1*length(t))
          var[[1]][[k]][curve,tp]<-var[[1]][[k]][curve,tp]+shift_val*rep(rgd,length(tp)) 
        }
      }
      #  else if (contamination_type=="weakly peak magnitude outlier")
      # {
      #   for (curve in outlier_nb)
      #   {
      #     shift_val<-ifelse(runif(1,-1,1)>0,1/4,-1/4)
      #     rv<-sample(1:length(t),1)
      #     st<-ifelse(rv>0.9*length(t),floor(0.9*length(t))-2,rv)
      #     tp<-st:(st+0.1*length(t))
      #     var[[k]][curve,tp]<-var[[k]][curve,tp]+shift_val*rep(rgd,length(tp)) 
      #   }        
      # 
       else if (contamination_type == "partial magnitude outlier")
      {
        for (curve in outlier_nb)
        {
          shift_val<-ifelse(runif(1,-1,1)>0,1/2,-1/2)
          st<-sample(2:((1-0.5)*length(t)-1),1)
          tp<-st:t[nbtime]
          var[[1]][[k]][curve,tp]<-var[[1]][[k]][curve,tp]+shift_val*rep(rgd,length(tp)) 
        }
      }
      # else if (contamination_type=="weakly partial magnitude outlier")
      # {
      #   for (curve in outlier_nb)
      #   {
      #     shift_val<-ifelse(runif(1,-1,1)>0,1/4,-1/4)
      #     st<-sample(2:((1-0.5)*length(t)-1),1)
      #     tp<-st:t[nbtime]
      #     var[[k]][curve,tp]<-var[[k]][curve,tp]+shift_val*rep(rgd,length(tp)) 
      #   }      
      # }
    }
  }
  invisible(var)
}
#######################################################
cluster_generation_minprop <- function(nbvar, n, nbtime, contamination_level, contamination_type, clusternumber, variat=c("amplitude","phase","shift"), partitionidentical)
{
  ## contamination_type:c("shape outlier I","shape outlier II",
  ## "shape outlier III","pure magnitude outlier","weakly pure magnitude outlier",
  ### "peak magnitude outlier","weakly peak magnitude outlier","partial magnitude outlier",
  ### "weakly partial magnitude outlier")
  outlier_nb<-sample(1:n,n*contamination_level)
  t<-seq(0,1,length.out=nbtime)

  var<-generate_cluster_minprop(nbvar,n,nbtime,clusternumber,variat,partitionidentical)
  var_outl<-introduce_outlier(var, outlier_nb, contamination_type)
  ###########################################  
  result<-list(var_outl[[1]],outlier_nb,label=var_outl[[2]])
  return (result)
}

############## plot the whole cluster
plot_clust <- function(result,rotate_angle, sc)
  ## result is the result of function cluster_generation
{
  var<-result[[1]]
  nbvar<-length(var)
  outlier_nb<-result[[2]]
  par(mfrow=c(1,ifelse(nbvar<=3,nbvar+1,3)), mar=c(5,4,4,1))
  for (k in 1:nbvar)
  {
    plot(t,var[[k]][1,],ylim=range(var[[k]]),type="n",ylab="Values",xlab="Time",main=paste("Variable ",k), cex.main=1.3, cex.lab = 1.3 , cex.axis = 1.3)
    
    lapply((1:n)[setdiff(1:n,outlier_nb)], function (ll) {
      for (cc in 1:clusternumber){
        if (cc<clusternumber& round((cc-1)*n/clusternumber)<ll &ll<=round((cc)*n/clusternumber))
        {
          lines(t,var[[k]][ll,],col=cc)
        } 
        else if (cc==clusternumber & (clusternumber-1)*round(n/clusternumber)<ll)
        {
          lines(t,var[[k]][ll,],col=cc)
        } 
        
      }
    })
    if (length(outlier_nb)>0)
    {
      apply(var[[k]][outlier_nb,],1, function (ll) {lines(t,ll,col="orange",lty=2)
      })
      text(t[1],var[[k]][outlier_nb,1],labels=outlier_nb,col="purple",cex=0.7)
    }
  }
  if (nbvar==3)
  {
    s3d<-scatterplot3d(x=var[[1]][1,],
                       y=var[[2]][1,],
                       z=var[[3]][1,],xlim=range(var[[1]]),ylim=range(var[[2]]),zlim=range(var[[3]])
                       ,type="n",xlab="Variable 1",ylab="Variable 2",zlab="Variable 3",
                       main="Trivarate Variables",  grid=TRUE,
                       scale.y=sc, angle=rotate_angle, cex.main=1.3, cex.lab = 0.9, cex.axis = 0.9)
    for (cc in 1:clusternumber)
    {
      if (cc<clusternumber)
      {
        normal_curve_index <- setdiff ((round((cc-1)*n/clusternumber)+1):round((cc)*n/clusternumber), outlier_nb)
        lapply(normal_curve_index, function(rowindex){
          s3d$points3d(var[[1]][rowindex,],var[[2]][rowindex,],var[[3]][rowindex,],type="l",col=cc)
        })
      } else
      {
        normal_curve_index <- setdiff(((clusternumber-1)*round(n/clusternumber)+1):n, outlier_nb)
        lapply(normal_curve_index,function(rowindex){
          s3d$points3d(var[[1]][rowindex,],var[[2]][rowindex,],var[[3]][rowindex,],type="l",col=cc)
        }) 
      }
    }
    if (length(outlier_nb)>0)
    {
      lapply(outlier_nb, function (ll) {s3d$points3d(var[[1]][ll,],var[[2]][ll,],var[[3]][ll,],col="orange",lty=2,type="l")
      })
      #text(var[[1]][outlier_nb,1],var[[2]][outlier_nb,1],var[[3]][outlier_nb,1],labels=outlier_nb,col="purple",cex=0.7)
    }
  }
  invisible(nbvar)
}

#### plot one representative from the cluster
plot_rep_line <- function(result,rotate_angle,sc)
  ## result is the result of function cluster_generation
{
  var<-result[[1]]
  nbvar<-length(var)
  n <- nrow(var[[1]])
  nbtime <- ncol(var[[1]])
  outlier_nb<-result[[2]]
  par(mfrow=c(1,ifelse(nbvar<=3,nbvar+1,3)), mar=c(5,4,4,1))
  for (k in 1:nbvar)
  {
    plot(t,var[[k]][1,],ylim=range(var[[k]]),type="n",ylab="Values",xlab="Time",main=paste("Variable ",k), cex.main=1.3, cex.lab = 1.3 , cex.axis = 1.3)
    
    lapply((1:n)[setdiff(1:n,outlier_nb)], function (ll) {
      for (cc in 1:clusternumber){
        if (cc<clusternumber&& round((cc-1)*n/clusternumber)<ll &&ll<=round((cc)*n/clusternumber))
        {
          lines(t,var[[k]][ll,],col=cc)
        } 
        else if (cc==clusternumber && (clusternumber-1)*round(n/clusternumber)<ll)
        {
          lines(t,var[[k]][ll,],col=cc)
        } 
        
      }
    })
    if (length(outlier_nb)>0)
    {
      apply(var[[k]][outlier_nb,],1, function (ll) {lines(t,ll,col="orange",lty=2)
      })
      text(t[1],var[[k]][outlier_nb,1],labels=outlier_nb,col="purple",cex=0.7)
    }
  }
 if (nbvar==3)
  {
    s3d<-scatterplot3d(x=var[[1]][1,],
                       y=var[[2]][1,],
                       z=var[[3]][1,],xlim=range(var[[1]]),ylim=range(var[[2]]),zlim=range(var[[3]])
                       ,type="n",xlab="Variable 1",ylab="Variable 2",zlab="Variable 3",
                       main="Trivarate Variables", grid=TRUE,
                        scale.y=sc, angle=rotate_angle , cex.main=1.3, cex.lab = 0.9, cex.axis = 0.9)
    for (cc in 1:clusternumber)
    {
      if (cc<clusternumber)
      {
        normal_curve_index <- setdiff ((round((cc-1)*n/clusternumber)+1):round((cc)*n/clusternumber), outlier_nb)
      } else
      {
        normal_curve_index <- setdiff(((clusternumber-1)*round(n/clusternumber)+1):n, outlier_nb)
      }
      rowindex<-normal_curve_index[1:2]
      s3d$points3d(var[[1]][rowindex,],var[[2]][rowindex,],var[[3]][rowindex,],type="l",col=cc)
      s3d$points3d(var[[1]][rowindex,1],var[[2]][rowindex,1],var[[3]][rowindex,1],type="p",pch =16, cex=1.2, col=cc)
      s3d$points3d(var[[1]][rowindex,nbtime],var[[2]][rowindex,nbtime],var[[3]][rowindex,nbtime],type="p", pch =17, cex=1.2, col=cc)
      
    }
    if (length(outlier_nb)>0)
    {
     ll <- outlier_nb[1]
     s3d$points3d(var[[1]][ll,],var[[2]][ll,],var[[3]][ll,],col="orange",lty=2,type="l")
     # text(var[[1]][ll,1],var[[2]][ll,1],var[[3]][ll,1],labels=ll,col="purple",cex=0.7)
    }
  }
  invisible(nbvar)
}
