function P_v2phi_nu,	nu, nu0, $
							ATOMICMASS = atomicmass, TEMP = temp, $
							; for maxwell distribution
							MEAN = mean, SD = sd
							; for gaussian distribution 

; convert maxwell-boltzmann velocity distribution to frequency distribution
; for a certain frequency of a given transition

  @fc                         ; load funcamental constants
  
  if keyword_set(TEMP) then begin
    m_atomic=atomicmass
    T=temp;
    phi = P_v_maxwell((nu0^2-nu^2)/(nu0^2+nu^2)*c, mu, sigma) $
    *(-4*nu0^2*nu)/(nu0^2+nu^2)^2*c
  endif
  if keyword_set(SD) then begin
    mu=mean
    sigma=sd
    phi = P_v_gaussian((nu0^2-nu^2)/(nu0^2+nu^2)*c, mu, sigma) $
    *(-4*nu0^2*nu)/(nu0^2+nu^2)^2*c
  endif
  return, phi

end
