function P_v_gaussian, v, mu, sigma
;+
;compute the probability DENSITY of a guassian distribution
;for a given velocity (normalized)
;
;v      : velocity
;mu     : mean/expectation value of veocity
;sigma  : standard deviation
;-
  @fc                 ; load funcamental constants
  
  P_v_gaussian = 1d/(sigma*sqrt(2d*!pi))*exp(-(1d/2d)*((v-mu)/sigma)^2d)
  return, P_v_gaussian

end
