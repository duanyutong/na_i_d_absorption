function lambda2v, lambda, lambda0
  @fc
  v=(lambda^2d - lambda0^2d)/(lambda^2d + lambda0^2d)*c
  return, v
end