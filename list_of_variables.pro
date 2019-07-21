  vdt = 1                             ; velocity distribution type: 1 for maxwellian, 2 for gaussian
  restrict_boundary_v = 0             ; choose whether or not impose a boundary condition for v
  restrict_boundary_density = 0       ; for density
  restrict_boundary_temp = 0          ; for temperature
  ; ----------- observables ----------
  r0 = 10d3*pc                        ; r_max. actual radius of the galaxy. measured horizontally
  r_proj = 0.5*r0                     ; projected radius. measured vertically
  inclination = acos(r_proj/r0)       ; calculate inclination angle. assume a perfectly circular disk
  ; ----------- dependences ----------
  ; for gaussian velocity distribution
  sigma = 150d3                       ; standard deviation
  mu = -1d3                           ; mean
  ; for maxwellian velocity distribution
  v_bulk_exp = 1                      ; [0, 3]
  number_density_exp = -1             ; [-3, 0]
  temp_exp = -1                       ; [-3, 0]
  v_max = 1d6               ; v_max at r_max; 
  v_bulk_coefficient = v_max/(r0)^v_bulk_exp ; automatically set and overrides manual v_bulk_coefficient
  number_density_max = 1d6  ;
  number_density_coefficient = number_density_max/(r0)^number_density_exp
  temp_max = 1d6
  temp_coefficient = temp_max/(r0)^temp_exp 
  ; -------- plotting settings ------- /without considering velocity due to hubble flow
  lambda_min = 5750d-10               ; lower limit of line profile domain
  lambda_max = 5950d-10               ; upper limit of line profile domain
  delta_lambda = 0.05d-10             ; increment of wavelength
  v_thermal_min = -10d3               ; lower limit of velocity distribution
  v_thermal_max = +10d3               ; upper limit of velocity distrubition
  delta_v_thermal = 0.01d3            ; increment of velocity 
  I0 = 1d                             ; original intensity of light beam 
  I_lambda0 = I0                      ; use a function instead for a radiation profile that depends on wavelength. eg. blackbody
  ; ----- transition 1 parameters ----
  f1 = 0.640511d                      ; oscillator strength
  bkbj1 = 1d                          ; the ratio b_k/b_j is 1 for bound-bound transition
  lambda1 = 5889.950954d-10           ; central wavelength for transition 1. 5890A in our case
  m_atomic1 = m_atomic_Na             ; atomic mass of atom 1
  ; ----- transition 2 parameters ----
  f2 = 0.319913d                      ; oscillator strength
  bkbj2 = 1d                          ; the ratio b_k/b_j is 1 for bound-bound transition
  lambda2 = 5895.924237d-10           ; central wavelength for transition 2. 5896A in our case
  m_atomic2 = m_atomic_Na             ; atomic mass of atom 2
  
  ; ===================================== main program ================================================
  ; ---- create wavelength array -----
  lambda_index = dindgen(ceil((lambda_max - lambda_min)/delta_lambda))
  lambda = [lambda_min + lambda_index*delta_lambda, lambda_max]
  lambda_n = size(lambda, /n_elements)
  ;------------------------------------------------------------------------------------
  ; set initial values for optical depth
  tau1_lambda = 0
  tau2_lambda = 0
  tau0_lambda = 0
    r_mean = r + delta_r/2d
    v_bulk = v_bulk(r_mean, v_bulk_coefficient, v_bulk_exp)
    T = model2_T(r_mean)
    ; compute integrated absoption coefficient
    s1 = acs_int(f1, bkbj1, lambda1, T)
    s2 = acs_int(f2, bkbj2, lambda2, T)
    case VDT of
      1:  begin       ; maxwell distribution
            v_thermal_index = dindgen(ceil((v_thermal_max - v_thermal_min)/delta_v_thermal))
            v_thermal = [v_thermal_min + v_thermal_index*delta_v_thermal, v_thermal_max]
            P_v1 = P_v_maxwell(v_thermal, m_atomic1, T)
            P_v2 = P_v_maxwell(v_thermal, m_atomic2, T)
            phi1 =  P_v2phi(v2lambda(lambda2v(lambda, lambda1) - v_bulk, lambda1), $
              lambda1, ATOMICMASS=m_atomic1, TEMP=T)
            phi2 =  P_v2phi(v2lambda(lambda2v(lambda, lambda2) - v_bulk, lambda2), $
              lambda2, ATOMICMASS=m_atomic2, TEMP=T)
          end
      2:  begin      ; gaussian distribution
            v_thermal_index = dindgen(ceil((v_thermal_max - v_thermal_min)/delta_v_thermal))
            v_thermal = [v_thermal_min + v_thermal_index*delta_v_thermal, v_thermal_max]
            P_v1 = P_v_gaussian(v_thermal, mu, sigma)
            P_v2 = P_v_gaussian(v_thermal, mu, sigma)
            phi1 =  P_v2phi(lambda, lambda1, MEAN=mu, SD=sigma)
            phi2 =  P_v2phi(lambda, lambda2, MEAN=mu, SD=sigma)
          end
    endcase
    v = v_add(v_thermal, v_bulk)
    v_n = size(v, /n_elements)
    ;-----------------------------------------------------------------------------------------------
    ; restrict range of lambda to avoid double peak due to negative velocities and such
    lambda_lowerlimit1 = min(v2lambda(v,lambda1))
    lambda_upperlimit1 = max(v2lambda(v,lambda1))
    print, 'limits of lambda as defined by the v_thermal+v_bulk for transition1 is' $
       + '[' + string(lambda_lowerlimit1) + ',' + string( lambda_upperlimit1) + ']'
    lambda_lowerlimit2 = min(v2lambda(v,lambda2))
    lambda_upperlimit2 = max(v2lambda(v,lambda2))
    print, 'limits of lambda as defined by the v_thermal+v_bulk for transition2 is' $
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
    column_density = qsimp('number_density2', r, r+delta_r, /double, eps=1d3)
    ; optical depth
    tau1_lambda = tau1_lambda + s1_lambda*column_density    ; optical depth for a specific wavelength, lambda
    tau2_lambda = tau2_lambda + s2_lambda*column_density    ; optical depth for a specific wavelength, lambda
    tau0_lambda = tau0_lambda + s0_lambda*column_density    ; optical depth for a specific wavelength, lambda
  endfor
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
  r_index = dindgen(ceil((r_max - r_min)/delta_r))
  r_array = [r_min + r_index*delta_r, r_max]
  r_n = size(r, /n_elements)
  if keyword_set(suppress_gfx) then begin
  endif else begin
;    w=0       ; window index
;    window, w++     ; velocity distribution of Na
    ps_start, file='/users/duan/Documents/research/model2/n^'$
      +STRTRIM(number_density2_exp, 2)+'_v^'+STRTRIM(v_bulk_exp,2)+'.eps', $
      /encap
      !p.multi=[0, 3, 2, 0, 0]
      margin = [2, 2, 1, 0]
        ; graph 1
        cgplot, r_array/pc, number_density2(r_array)/1d6,$
          background='white',color='red',axiscolor='black',$
          title = 'Number Density of Na', $
          xtitle = 'Radius/pc', ytitle = 'Number density/'+textoidl('cm^{-3}')
        ; graph 2
        cgplot, r_array/pc, v_bulk(r_array, v_bulk_coefficient, v_bulk_exp)/1d3,$
          background='white', color='red', axiscolor='black', $
          title = 'Bulk Velocity', $
          xtitle = 'Radius/pc', ytitle = 'Velocity/'+textoidl('km\cdot')+textoidl('s^{-1}')
        ; graph 3
        cgplot, r_array/pc, model2_T(r_array), $
          background='white',color='red',axiscolor='black', $
          title = 'Temperature', $
          xtitle = 'Radius/pc', ytitle = 'Temperature/K'
        ; graph 4
        cgplot, v_thermal/1d3, P_v_maxwell(v_thermal, m_atomic1, mean(T)), $
          background='white', color='red',axiscolor='black', $
          title = 'Maxwellian Distribution (Thermal)', $
          xtitle = 'Velocity/'+textoidl('km\cdot')+textoidl('s^{-1}'), ytitle = 'Probability density'
        ; graph 5
        cgplot, lambda*1d10, tau1_lambda, /ylog, $
          background='white',color='red',axiscolor='black', charthick = 1, thick = 2, $
          title = 'Optical depth', $
          xtitle = 'Wavelength/'+cgsymbol('Angstrom'), ytitle = 'Optical depth', $
          xr = [5850d, 5925d], xsty=1
        cgoplot, lambda*1d10, tau2_lambda, /ylog, $
          background='white',color='blue',axiscolor='black', charthick = 1, thick = 2
        ; graph 6
        cgplot, lambda*1d10, I1_lambda/I_lambda0, $
          background='white',color='red',axiscolor='black', charthick = 1, thick = 2, $
          title = 'Relative Intensity', $
          xtitle = 'Wavelength/'+cgsymbol('Angstrom'), ytitle = 'Relative intensity', $
          xr = [5850d, 5925d], xsty=1
        cgoplot, lambda*1d10, I2_lambda/I_lambda0, $
          background='white',color='blue',axiscolor='black', charthick = 1, thick = 2
        cgoplot, lambda*1d10, I0_lambda/I_lambda0,$
          background='white',color='black',axiscolor='black', charthick = 1, thick = 2
      ; for transitions of two types of atoms with different veolocity distributions
      ; cgplot, v, P_v2, background='white', color='red', axiscolor='black'
      !p.multi=0
    ps_end
  endelse
  ;-------------------------------------------------------------------------------------
  ; save cache; debug
  window, 2
        cgplot, lambda*1d10, I1_lambda/I_lambda0, background='white',color='red',axiscolor='black', $
        title = 'Relative intensity'
      cgoplot, lambda*1d10, I2_lambda/I_lambda0, background='white',color='blue',axiscolor='black'
      cgoplot, lambda*1d10, I0_lambda/I_lambda0, background='white',color='black',axiscolor='black'
  save, /variables, filename='var.sav'
  stop
end
