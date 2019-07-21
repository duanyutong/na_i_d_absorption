function number_density2, r
  ; number of absorbing atoms in unit volume along the radial direction
  ; rho - r^-2
  common ndshare, number_density2_exp, number_density2_coefficient
  number_density = number_density2_coefficient*r^(number_density2_exp)
  return, number_density
end