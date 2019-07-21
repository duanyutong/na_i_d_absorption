function v2lambda, v, lambda0
  @fc
  beta = v/c
  lambda=lambda0*sqrt((1+beta)/(1-beta))
  return, lambda
end