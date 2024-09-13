
ft = function(time){
  time
}

nearest = function(x,u){
  which.min(abs(x-u))
}

# For constant rates, without ceiling
visitingPattern = function(rate,maxTime){
  times=NULL
  time = 0
  while(time<=maxTime){
    interval = rexp(1,rate=rate)
    
    time = time+interval
    if(time>=maxTime) break
    if(length(times)>500) {
      warnings("more than 500 visits.")
    }
    times = c(times, time)
  }
  return(times)
}

# For constant rates, with ceiling
visitingPattern2 = function(rate,maxTime){
  times=NULL
  time = 0
  while(time<=maxTime){
    interval = ceiling(rexp(1,rate=rate))
    
    time = time+interval
    if(time>=maxTime) break
    if(length(times)>500) {
      warnings("more than 500 visits.")
    }
    times = c(times, time)
  }
  return(times)
}

# For varying rates, with ceiling
visitingPattern3 = function(rate,maxTime){
  times=NULL
  time = 0
  i=1
  while(time<=maxTime){
    interval = ceiling(rexp(1,rate=rate[i])) # round up for multiple biomarkers
    time = time+interval
    i=time+1 # since data starts from time=0
    if(time>=maxTime) break
    if(length(times)>500) {
      warnings("more than 500 visits.")
    }
    times = c(times, time)
  }
  return(times)
}
