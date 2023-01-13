OptionPrice = function(S, K, sigma, r, q, ttm, type)
{
  b = r - q
  t = ttm/365
  
  d1 = (log(S / K) + (b + sigma ^ 2 / 2) * t) / (sigma * sqrt(t))
  d2 = d1 - sigma * sqrt(t)
  if(type == "call")
  {
    price = S * exp((b - r) * t) * pnorm(d1) - K * exp(-r * t) * pnorm(d2)
  }
  else if (type == "put")
  {
    price = (K * exp(-r * t) * pnorm(-d2) - S * exp((b - r) * t) * pnorm(-d1))
  }
  
  return(price)
}

OptionVega = function (S, K, sigma, r, q, ttm, type)
{
  b = r - q
  t = ttm/365
  d1 = (log(S / K) + (b + sigma ^ 2 / 2) * t) / (sigma * sqrt(t))
  return  (S * sqrt(t) * dnorm(d1))
}

OptionIV = function (S, K, price, r, q, ttm, type)
{
  maxIter = 100
  precision = 0.000000001
  
  vol = 0.20 #initial guess
  for (i in 1:maxIter)
  {
    value = OptionPrice (S=S, K=K, sigma=vol, r=r, q=q, ttm=ttm, type=type)
    vega = OptionVega (S=S, K=K, sigma=vol, r=r, q=q, ttm=ttm, type=type)
    
    if (vega==0)
      return (0.00001)
    
    diff = price - value

    if (abs(diff) < precision)
      return (vol)
    else 
      vol = vol + 0.5 * diff / vega
   }
  
  return(vol)
}

OptionGamma = function (S, K, sigma, r, q, ttm, type)
{
  b = r - q
  t = ttm/365
  
  d1 = (log(S / K) + (b + sigma ^ 2 / 2) * t) / (sigma * sqrt(t))
  return (dnorm(d1) / (S * sigma * sqrt(t)))
}

OptionVolga = function(S, K, sigma, r, q, ttm, type)
{
  b = r - q
  t = ttm/365
  
  d1 = (log(S / K) + (b + sigma ^ 2 / 2) * t) / (sigma * sqrt(t))
  d2 = d1 - sigma * sqrt(t)
  
  return (S * sqrt(t) * dnorm(d1) * d1 * d2 / sigma)
}

OptionRho = function(S, K, sigma, r, q, ttm, type)
{
  b = r - q
  t = ttm/365
  
  d1 = (log(S / K) + (b + sigma ^ 2 / 2) * t) / (sigma * sqrt(t))
  d2 = d1 - sigma * sqrt(t)
  
  if(type == "call")
    rho =  K * t * exp(-r*t) * pnorm (d2)
  else if (type == "put")
    rho = -K * t * exp(-r*t) * pnorm (-d2)
  
  return (rho)
}

OptionDelta = function (S, K, sigma, r, q, ttm, type)
{
  b = r - q
  t = ttm/365
  
  d1 = (log(S / K) + (b + sigma ^ 2 / 2) * t) / (sigma * sqrt(t))
  
  if(type == "call")
    delta = exp(-q * t) * pnorm(d1)
  else if (type == "put")
    delta =  exp(-q * t) * (pnorm(d1)-1)
  
  return(delta)
}

OptionTheta = function (S, K, sigma, r, q, ttm, type)
{
  b = r - q
  t = ttm/365
  
  d1 = (log(S / K) + (b + sigma ^ 2 / 2) * t) / (sigma * sqrt(t))
  d2 = d1 - sigma * sqrt(t)
  
  if(type == "call")
    theta = -S*dnorm(d1)*sigma*exp(-q*t) / (2 * sqrt(t)) + q*S*pnorm(d1)*exp(-q*t) - r*K*exp(-r*t)*pnorm(d2)
  else if (type == "put")
    theta = -S*dnorm(d1)*sigma*exp(-q*t) / (2 * sqrt(t)) - q*S*pnorm(d1)*exp(-q*t) + r*K*exp(-r*t)*pnorm(d2)

  return(theta)   
}