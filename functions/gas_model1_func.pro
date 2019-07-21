function gas_model1_func, density_coefficient, SUPPRESS_GFX=suppress_gfx
  ;+
  ; function version of simple model for a parcel of gas
  ;-
  @fc                          ; load funcamental constants
  ; define input variables
  ; velocity distribution type. 1 for maxwellian, 2 for gaussian
  ; var for maxwellian dist: m_atomic, T
  ; var for guassian dist: mu, sigma
  VDT = 2
  mu = 1d                            ; mean of gaussian velocity distribution
  sigma = 150d3                     ; standard deviation of gaussian velocity distribution
  lambda_min = 5870d-10                ; lower limit of line profile domain
  lambda_max = 5920d-10                 ; upper limit of line profile domain
  delta_lambda = 0.01d-10               ; increment of wavelength
  v_min = -1000d3                      ; lower limit of velocity distribution
  v_max = +1000d3                     ; upper limit of velocity distrubition
  delta_v = 1d3                   ; increment of velocity
  T = 30d6                          ; temperature of the gas, to be used for calculating acs_int & p_v_maxwell
  I0 = 1d                         ; original intensity of light beam 
  I_lambda0 = I0                ; use a function instead for a radiation profile that depends on wavelength. eg. blackbody
  r_min = 0d                      ; lower limit of line of sight
  r_max = 1d-7*pc                       ; upper limit of line of sight
; for transition 1--------------------------------------------------------------------
  f1 = 0.640511d                         ; oscillator strength
  bkbj1 = 1d                      ; the ratio b_k/b_j is 1 for bound-bound transition
  lambda1 = 5889.950954d-10                    ; central wavelength for transition 1. 5890A in our case
  m_atomic1 = m_atomic_Na                  ; atomic mass of atom 1
  ; for transition 2--------------------------------------------------------------------
  f2 = 0.319913d                         ; oscillator strength
  bkbj2 = 1d                      ; the ratio b_k/b_j is 1 for bound-bound transition
  lambda2 = 5895.924237d-10                    ; central wavelength for transition 2. 5896A in our case
  m_atomic2 = m_atomic_Na                  ; atomic mass of atom 2
  ;------------------------------------------------------------------------------------
  lambda_index = dindgen(ceil((lambda_max - lambda_min)/delta_lambda))
  lambda = [lambda_min + lambda_index*delta_lambda, lambda_max]
  lambda_n = size(lambda, /n_elements)
  ;------------------------------------------------------------------------------------
  ; compute integrated absoption coefficient
  s1 = acs_int(f1, bkbj1, lambda1, T)
  s2 = acs_int(f2, bkbj2, lambda2, T)
  ;------------------------------------------------------------------------------------
  ; wavelength distribution
  case VDT of
    ; maxwell distribution
    1:  begin
          v_index = dindgen(ceil((v_max - v_min)/delta_v))
          v = [v_min + v_index*delta_v, v_max]
          v_n = size(v, /n_elements)
          P_v1 = P_v_maxwell(v, m_atomic1, T)
          P_v2 = P_v_maxwell(v, m_atomic2, T)
          phi1 =  P_v2phi(lambda, lambda1, ATOMICMASS=m_atomic1, TEMP=T)
          phi2 =  P_v2phi(lambda, lambda2, ATOMICMASS=m_atomic2, TEMP=T)
        end
    ; gaussian distribution
    2:  begin
          v_index = dindgen(ceil((v_max - v_min)/delta_v))
          v = [v_min + v_index*delta_v, v_max]
          v_n = size(v, /n_elements)
          P_v1 = P_v_gaussian(v, mu, sigma)
          P_v2 = P_v_gaussian(v, mu, sigma)
          phi1 =  P_v2phi(lambda, lambda1, MEAN=mu, SD=sigma)
          phi2 =  P_v2phi(lambda, lambda2, MEAN=mu, SD=sigma)
        end
  endcase
  ;-----------------------------------------------------------------------------------------------
  ; restrict range of lambda to avoid double peak due to negative velocities and such
  lambda_lowerlimit1 = min(v2lambda(v,lambda1))
  lambda_upperlimit1 = max(v2lambda(v,lambda1))
  print, 'limits of lambda as defined by the velocity distribution for transition1 is' $
     + '[' + string(lambda_lowerlimit1) + ',' + string( lambda_upperlimit1) + ']'
  lambda_lowerlimit2 = min(v2lambda(v,lambda2))
  lambda_upperlimit2 = max(v2lambda(v,lambda2))
  print, 'limits of lambda as defined by the velocity distribution for transition2 is' $
     + '[' + string(lambda_lowerlimit2) + ',' + string( lambda_upperlimit2) + ']'
  lambda_lowerlimit1_index = value_locate(lambda, lambda_lowerlimit1)
  lambda_upperlimit1_index = value_locate(lambda, lambda_upperlimit1)
  lambda_lowerlimit2_index = value_locate(lambda, lambda_lowerlimit2)
  lambda_upperlimit2_index = value_locate(lambda, lambda_upperlimit2)
  if lambda_lowerlimit1 lt lambda_min then begin
    print, 'decrease lambda_min to display all data. execution halted'
    stop      ; lambda_lowerlimit outside the lambda interval to be plotted on x-axis
  endif
  if lambda_upperlimit1 gt lambda_max then begin
    print, 'increase lambda_max to display all data. execution halted'
    stop      ; lambda_upperlimit on the right of the interval
  endif
  if lambda_lowerlimit2 lt lambda_min then begin
    print, 'decrease lambda_min to display all data. execution halted'
    stop      ; lambda_lowerlimit outside the lambda interval to be plotted on x-axis
  endif
  if lambda_upperlimit2 gt lambda_max then begin
    print, 'increase lambda_max to display all data. execution halted'
    stop      ; lambda_upperlimit on the right of the interval
  endif
  phi1(0:lambda_lowerlimit1_index)=0       ; null meaningless data outside what's defined by the velocity distribution
  phi1(lambda_upperlimit1_index: size(phi1, /n_elements)-1)=0
  phi2(0:lambda_lowerlimit2_index)=0       ; null meaningless data outside what's defined by the velocity distribution
  phi2(lambda_upperlimit2_index: size(phi2, /n_elements)-1)=0
  ;---------------------------------------------------------------------------------------
  ; absorption cross section for a specific wavelength, lambda
  s1_lambda = s1*phi1
  s2_lambda = s2*phi2
  s0_lambda = s1_lambda + s2_lambda
  ; column density
  column_density = density_coefficient*qsimp('number_density', r_min, r_max, /double)
  ; optical depth
  tau1_lambda = s1_lambda*column_density    ; optical depth for a specific wavelength, lambda
  tau2_lambda = s2_lambda*column_density    ; optical depth for a specific wavelength, lambda
  tau0_lambda = s0_lambda*column_density    ; optical depth for a specific wavelength, lambda
  ; resulting intensity at a specific wavelength
  I1_lambda = I_lambda0*exp(-tau1_lambda)
  I2_lambda = I_lambda0*exp(-tau2_lambda)
  I0_lambda = I_lambda0*exp(-tau0_lambda)
  ;---------------------------------------------------------------------------------------------------
  ; equivalent width
  equivalent_width1 = discrete_area(lambda, 1-I1_lambda/I_lambda0)
  equivalent_width2 = discrete_area(lambda, 1-I2_lambda/I_lambda0)
  equivalent_width0 = discrete_area(lambda, 1-I0_lambda/I_lambda0)
  ; plots---------------------------------------------------------------------------------
  if keyword_set(suppress_gfx) then begin
  endif else begin
    w=0       ; window index
    window, w++     ; velocity distribution of Na
    cgplot, v, P_v1, background='white',color='red',axiscolor='black'
    ; for transitions of two types of atoms with different veolocity distributions
    ; cgplot, v, P_v2, background='white', color='red', axiscolor='black'
    window, w++     ; optical depth
    cgplot, lambda, phi1, background='white',color='red',axiscolor='black'
    window, w++             ; absorption line profile for transition1
    cgplot, lambda, I1_lambda/I_lambda0, background='white',color='red',axiscolor='black'
    window, w++             ; absorption line profile for transition2
    cgplot, lambda, I2_lambda/I_lambda0, background='white',color='red',axiscolor='black'
    window, w++             ; absorption line profile for doublet
    cgplot, lambda, I0_lambda/I_lambda0, background='white',color='red',axiscolor='black'
  endelse
  ;-------------------------------------------------------------------------------------
  ; save cache; debug
  save, /variables, filename='var.sav'
  return, [column_density, equivalent_width1, equivalent_width2]
  
end