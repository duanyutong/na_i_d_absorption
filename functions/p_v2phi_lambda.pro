function P_v2phi_lambda, lambda, lambda0, $
                    ATOMICMASS = atomicmass, TEMP = temp, $
					; for maxwell distribution
                    MEAN = mean, SD = sd
					; for gaussian distribution 

;convert maxwell-boltzmann velocity distribution to wavenlength distribution
;for a certain trasition and given wavelength to be observed

  @fc                         ; load funcamental constants
  
  if keyword_set(TEMP) then begin
    m_atomic=atomicmass
    T=temp;
    phi = $
		P_v_maxwell((lambda^2-lambda0^2)/(lambda^2+lambda0^2)*c, m_atomic, T)$
		*(4*lambda0^2*lambda)/(lambda0^2+lambda^2)^2*c
  endif
  if keyword_set(SD) then begin
    mu=mean
    sigma=sd
    phi = $
		P_v_gaussian((lambda^2-lambda0^2)/(lambda^2+lambda0^2)*c, mu, sigma)$
		*(4*lambda0^2*lambda)/(lambda0^2+lambda^2)^2*c
  endif  
  return, phi

end
