function P_v_maxwell, v, m_atomic, T
; maxwell-blotzmann velocity distrition in one dimension
; compute the probability DENSITY for a given velocity (normalized)
;
; m_atomic : atomic mass
;           m=22.989769282 for sodium
; T        : temperature of gas
; v        : velocity

  @fc                         ; load funcamental constants
  
  m=m_u*m_atomic              ;compute actual mass of the atom
  P_v_maxwell = sqrt(m/(2*!pi*k*T))*exp(-(m*v^2)/(2*k*T))
  return, P_v_maxwell

end
