function gas_model3_func, n, SUPPRESS_GFX=suppress_gfx
  @fc                                 ; load funcamental constants
  common shared_parameters, number_density_coefficient, number_density_exp, $
                          temp_coefficient,           temp_exp, $
                          v_bulk_coefficient,         v_bulk_exp          ; share parameters with functions
  ; =================================== define input variables ===========================================
  ; variables for maxwellian distribution : m_atomic, T
  ; variables for guassian distribution   : mu, sigma
  ; Geometry assumptions:   observing along z direction
  ;                         major axis of the ellipse observed always oerlaps the x-axis
  ;                         galactic center is at (0,0) on x-y plane
  ; ----------- model options --------
  number_density_coefficient = n
  spot = [0, 0]                       ; location of the pixel (grid unit) observed in x-y plane
  grid_radius = 10d                   ; number of pixels along radius. half length of each dimension
  vdt = 2                             ; velocity distribution type: 1 for maxwellian , 2 for gaussian
  restrict_boundary_v_bulk = 0        ; choose whether or not impose a boundary condition for v
  restrict_boundary_density = 0       ; for density
  restrict_boundary_temp = 0          ; for temperature
  ; ----------- observables ----------
  ; we can measure original/projected area or radius
  r0 = 10d3*pc                        ; r_max. actual radius of the galaxy. measured horizontally
  r_proj = 0.5*r0                     ; projected radius. measured vertically
  inclination = 2                     ; inclination angle between galactic disc and z+ axis. range [0, pi)
  v_hubble = 0
  ; ----------- dependences ----------
  sigma = 100d3                       ; standard deviation   ; for gaussian velocity distribution
  mu = 0d                             ; mean
  number_density_exp = -1             ; [-3, 0]
  number_density_max = 1d6            ; max density at r_max; 
  temp_exp = -1                       ; [-3, 0]
  temp_max = 1d5                      ; for maxwell velocity distribution
  v_bulk_exp = 1                      ; [ 0, 3]
  v_bulk_max = 1d6
  ; -------- plotting settings ------- /without considering velocity due to hubble flow
  lambda_min = 5800d-10               ; lower limit of line profile domain in the plot and calculation
  lambda_max = 5950d-10               ; upper limit of line profile domain
  delta_lambda = 0.10d-10             ; increment of wavelength
  v_gaussian_max = 800d3              ; upper limit of thermal speed
  delta_v_gaussian = 0.05d3           ; increment of thermal velocity 
  I0 = 1d                             ; original intensity of light beam 
  I_lambda0 = I0                      ; use a function instead for a radiation profile that depends on wavelength. eg. blackbody
  ; ----- Na I D1 transition parameters ----
  f1 = 0.640511d                      ; oscillator strength
  bkbj1 = 1d                          ; the ratio b_k/b_j is 1 for bound-bound transition
  lambda1 = 5897.56661715d-10         ; central wavelength for transition 1. 5890A in our case
  m_atomic1 = m_atomic_Na             ; atomic mass of atom 1
  ; ----- Na I D2 ----
  f2 = 0.319913d                      ; oscillator strength
  bkbj2 = 1d                          ; the ratio b_k/b_j is 1 for bound-bound transition
  lambda2 = 5891.58326415d-10         ; central wavelength for transition 2. 5896A in our case
  m_atomic2 = m_atomic_Na             ; atomic mass of atom 2
  
  ;=========================== set coefficients of dependence relations  ==========================
  case restrict_boundary_temp of
      0:  begin
            temp_coefficient = 1d18
          end
      1:  begin
            temp_coefficient = temp_max/(r0)^temp_exp
          end
  endcase
  case restrict_boundary_v_bulk of
      0:  begin
            v_bulk_coefficient = 1d-13
          end
      1:  begin
            v_bulk_coefficient = v_bulk_max/(r0)^v_bulk_exp
          end
  endcase

  ; ===================================== construct grid ================================================
  ; allocate memory
  side = 2*grid_radius+1          ; total number of pixels along each dimension
  xhat = [1, 0, 0]
  yhat = [0, 1, 0]
  zhat = [0, 0, 1]
  normal = [0, sin(inclination-!pi/2), cos(inclination-!pi/2)]    ; disk's normal vector always in the y-z plane
  occupy = dblarr(side,side,side)     ; occupy cartesian space! we physicists are the 99%!
  grid = {x:              occupy, $   ; index of x coordinate
          y:              occupy, $   ; index of y coordinate
          z:              occupy, $   ; index of z coordinate
          r:              occupy, $   ; radial distance to galactic center
          theta:          occupy, $   ; angle to galactic disk plane
          number_density: occupy, $
          temp:           occupy, $
          v_bulk:         occupy, $
          v_bulk_x:       occupy, $
          v_bulk_y:       occupy, $
          v_bulk_z:       occupy, $
          plane:          occupy  $   ; use as logical field. if the plane lies in the block, value is 1
          }
  origin = [grid_radius+1,grid_radius+1,grid_radius+1]
  for i = 0, side-1 do begin         ; assign index numbers for x, y, z
    occupy = replicate(i-grid_radius,side,side)
    grid.x(i,*,*) = occupy
    grid.y(*,i,*) = occupy
    grid.z(*,*,i) = occupy
  endfor
  for i = 0, side-1 do begin
    for j = 0, side-1 do begin
      for k = 0, side-1 do begin         ; assume azimuthally symmetric wind
        position = [grid.x(i,j,k), grid.y(i,j,k), grid.z(i,j,k)]    ; dimensionless, simply index numbers for cartesian grid
        grid.r(i,j,k) = norm(position)*r0/grid_radius               ; calculate radial distance in unit of length (m)
        grid.theta(i,j,k) = acos(                                 $ 
                                sum(position*normal)              $ ; dot product
                                /(norm(position)*norm(normal))    $ ; product of norms
                                )
        position_spherical = [grid.r(i,j,k), grid.theta(i,j,k)]
        ; --------- fill in physical conditions ----------
        if grid.r(i,j,k) gt r0 then begin
          grid.number_density(i,j,k) = 0d         ; outside the spherical wind
        endif else begin
          grid.number_density(i,j,k) = number_density_spherical(position_spherical)
          grid.temp(i,j,k) = temp_spherical(position_spherical)
          v_bulk = position/norm(position)*abs(v_bulk_spherical(position_spherical))    
          ; v_bulk_spherical returns +- velocities due to spherical coordinate system cos(theta)
          grid.v_bulk(i,j,k) = norm(v_bulk)
          grid.v_bulk_x(i,j,k) = v_bulk(0)
          grid.v_bulk_y(i,j,k) = v_bulk(1)
          grid.v_bulk_z(i,j,k) = v_bulk(2)*(-1d)
          if keyword_set(suppress_gfx) then begin
          endif else begin
            print, 'finished filling grid(' + string(i)+','+string(j)+','+string(k)+')'
          endelse
        endelse
      endfor
      ; --------- determine if block is a plane block -----
      los_logic = grid.theta(i,j,*) gt !pi/2
      plane_index = SEARCH_ARRAY(los_logic, 0)
      plane_index = plane_index(0)
      grid.plane(i,j,plane_index) = 1
    endfor
  endfor
  ; --------- fill galactic center by perturbation --------------
  position = [0,1*pc,0]    ; dimensionless, simply index numbers for cartesian grid
  grid.r(grid_radius,grid_radius,grid_radius) = 1*pc
  grid.theta(grid_radius,grid_radius,grid_radius) = 1d-3
  position_spherical = [1*pc, 1d-3]
  grid.number_density(grid_radius,grid_radius,grid_radius) = number_density_spherical(position_spherical)
  grid.temp(grid_radius,grid_radius,grid_radius) = temp_spherical(position_spherical)
  v_bulk = position/norm(position)*abs(v_bulk_spherical(position_spherical))    
  ; v_bulk_spherical returns +- velocities due to spherical coordinate system cos(theta)
  grid.v_bulk(grid_radius,grid_radius,grid_radius) = norm(v_bulk)
  grid.v_bulk_x(grid_radius,grid_radius,grid_radius) = v_bulk(0)
  grid.v_bulk_y(grid_radius,grid_radius,grid_radius) = v_bulk(1)
  grid.v_bulk_z(grid_radius,grid_radius,grid_radius) = v_bulk(2)*(-1d)
  if keyword_set(suppress_gfx) then begin
  endif else begin
    print, 'finished filling grid(' + string(i)+','+string(j)+','+string(k)+')'
  endelse
  ; ===================================== calculate spaxels =============================================
  ; ------- create empty matrices --------------
  lambda_index = dindgen(ceil((lambda_max - lambda_min)/delta_lambda))
  lambda = [lambda_min + lambda_index*delta_lambda, lambda_max]     ; create wavelength array
  lambda_n = size(lambda, /n_elements)
  spaxel = dblarr(side,side,3,lambda_n)         ; first two dimensions for pixel position, forth dimension for storing spectrum
  tau = dblarr(side,side,3,lambda_n)            ; third dimsnsion for specifying tau0, tau1, and tau2
  column_density = dblarr(side,side)
  v_gaussian_index = dindgen(ceil((2*v_gaussian_max/delta_v_gaussian)))
  v_gaussian = [-v_gaussian_max + v_gaussian_index*delta_v_gaussian, v_gaussian_max]
  ; --------- integrate along line of sight ----------- 
;  for i = 0, side-1 do begin
;    for j = 0, side-1 do begin
      i = spot(0)+grid_radius
      j = spot(1)+grid_radius
      for k = side-1, 0, -1 do begin      ; integrate toward -z direction
        if (grid.plane(i,j,k) eq 1) then break
        ; -------compute integrated absoption coefficient ------------
        s1 = acs_int(f1, bkbj1, lambda1, grid.temp(i,j,k))
        s2 = acs_int(f2, bkbj2, lambda2, grid.temp(i,j,k))
        ; ---------- determine wavelengths of interest ------------
        v = v_add(                              $
                  grid.v_bulk_z(i,j,k),         $
                  v_add(v_gaussian, v_hubble)   $
                  )
        ; v_n = size(v, /n_elements)
        info_lambda1_min = min(v2lambda(v,lambda1))
        info_lambda2_min = min(v2lambda(v,lambda2))
        info_lambda1_max = max(v2lambda(v,lambda1))
        info_lambda2_max = max(v2lambda(v,lambda2))
        if keyword_set(suppress_gfx) then begin
        endif else begin
          print, 'wavelength range corresponding to defined gaussian, v_hubble and v_bulk_z for Na I D1 is'
          print, '[' + string(info_lambda1_min*1d10) + ',' + string( info_lambda1_max*1d10) + ']'
          print, 'wavelength range corresponding to defined gaussian, v_hubble and v_bulk_z for Na I D2 is'
          print, '[' + string(info_lambda2_min*1d10) + ',' + string( info_lambda2_max*1d10) + ']'
        endelse
        ; ---------- velocity and wavelength distribution ------------
        v_rest1 = v_add(                              $ ; rest gaussian velocity in the grid block frame
                      lambda2v(lambda, lambda1),      $
                      v_add(                          $
                            -grid.v_bulk_z(i,j,k),    $
                            -v_hubble                 $
                            )                         $
                      )
        v_rest2 = v_add(                              $ ; rest gaussian velocity in the grid block frame
                      lambda2v(lambda, lambda2),      $
                      v_add(                          $
                            -grid.v_bulk_z(i,j,k),    $
                            -v_hubble                 $
                            )                         $
                      )
        case vdt of
          1:  begin       ; maxwell distribution
                phi1 = P_v2phi(                               $
                              v2lambda(v_rest1, lambda1),     $
                              lambda1,                        $
                              ATOMICMASS=m_atomic1, TEMP=T    $
                              )                               
                phi2 = P_v2phi(                               $
                              v2lambda(v_rest2, lambda2),     $
                              lambda2,                        $
                              ATOMICMASS=m_atomic2, TEMP=T    $
                              )
              end
          2:  begin      ; gaussian distribution
                phi1 = P_v2phi(                               $
                              v2lambda(v_rest1, lambda1),     $
                              lambda1,                        $
                              MEAN=mu, SD=sigma               $
                              )   ; assume same gaussian for all locations
                phi2 = P_v2phi(                               $
                              v2lambda(v_rest2, lambda2),     $
                              lambda2,                        $
                              MEAN=mu, SD=sigma               $
                              )
              end
        endcase
        ; --------------------- absorption cross section for each specific wavelength ---------------------
        s1_lambda = s1*phi1
        s2_lambda = s2*phi2
        s0_lambda = s1_lambda + s2_lambda
        ; column density
        column_density_block = grid.number_density(i,j,k)*r0/grid_radius
        ; optical depth
        tau(i,j,0,*) = tau(i,j,0,*) + s0_lambda*column_density_block
        tau(i,j,1,*) = tau(i,j,1,*) + s1_lambda*column_density_block
        tau(i,j,2,*) = tau(i,j,2,*) + s2_lambda*column_density_block
        column_density(i,j) = column_density(i,j) + column_density_block
        if keyword_set(suppress_gfx) then begin
        endif else begin
          print, 'integrated block z='+ string(grid.z(i,j,k))
        endelse
      endfor
;    endfor
;  endfor
  ; ---------------- resulting intensity ---------------- 
  spaxel(i,j,0,*) = I_lambda0*exp(-tau(i,j,0,*))
  spaxel(i,j,1,*) = I_lambda0*exp(-tau(i,j,1,*))
  spaxel(i,j,2,*) = I_lambda0*exp(-tau(i,j,2,*))
  ; ---------------- equivalent width ------------------
  equivalent_width0 = discrete_area(lambda, 1-spaxel(i,j,0,*)/I_lambda0)
  equivalent_width1 = discrete_area(lambda, 1-spaxel(i,j,1,*)/I_lambda0)
  equivalent_width2 = discrete_area(lambda, 1-spaxel(i,j,2,*)/I_lambda0)
  ; ---------------- plots ----------------------------
  if keyword_set(suppress_gfx) then begin
  endif else begin
;  r_index = dindgen(ceil((r_max - r_min)/delta_r))
;  r_array = [r_min + r_index*delta_r, r_max]
;  r_n = size(r, /n_elements)
;    w=0       ; window index
;    window, w++     ; velocity distribution of Na
    z = grid.z(i,j,*)*r0/grid_radius
    ps_start, file='/users/duan/Documents/research/model3/n^' $
                    +strtrim(number_density_exp,2)+'_v^' $
                    +strtrim(v_bulk_exp,2)+'_T^' $
                    +strtrim(temp_exp,2)+'.eps', /encap
      !p.multi=[0, 3, 2, 0, 0]
      margin = [2, 2, 1, 0]
      ; graph 1
      cgplot, z/pc, grid.number_density(i,j,*)/1d6, /ylog, $
        background='white',color='red',axiscolor='black',$
        title = 'Number Density', $
        xtitle = 'Line of Sight/pc', ytitle = 'Number density/'+textoidl('cm^{-3}')
      ; graph 2
      cgplot, z/pc, grid.v_bulk_z(i,j,*)/1d3,$
        background='white', color='red', axiscolor='black', $
        title = 'Bulk Velocity', $
        xtitle = 'Line of Sight/pc', ytitle = 'Velocity/'+textoidl('km\cdot')+textoidl('s^{-1}')
      ; graph 3
      cgplot, z/pc, grid.temp(i,j,*), /ylog, $
        background='white',color='red',axiscolor='black', $
        title = 'Temperature', $
        xtitle = 'Line of Sight/pc', ytitle = 'Temperature/K'
      ; graph 4
      cgplot, v_gaussian/1d3, P_v_gaussian(v_gaussian, mu, sigma), $
        background='white', color='red',axiscolor='black', $
        title = 'Gaussian Distribution', $
        xtitle = 'Velocity/'+textoidl('km\cdot')+textoidl('s^{-1}'), ytitle = 'Probability density'
      ; graph 5
      cgplot, lambda*1d10, tau(i,j,1,*), /ylog, $
        background='white',color='red',axiscolor='black', charthick = 1, thick = 2, $
        title = 'Optical depth', $
        xtitle = 'Wavelength/'+cgsymbol('Angstrom'), ytitle = 'Optical depth', $
        xr = [5850d, 5910d], xsty=1
      cgoplot, lambda*1d10, tau(i,j,2,*), /ylog, $
        background='white',color='blue',axiscolor='black', charthick = 1, thick = 2
      ; graph 6
      cgplot, lambda*1d10, spaxel(i,j,0,*)/I_lambda0, $
        background='white',color='black',axiscolor='black', charthick = 1, thick = 2,  $
        title = 'Relative Intensity', $
        xtitle = 'Wavelength/'+cgsymbol('Angstrom'), ytitle = 'Relative intensity', $
        xr = [5850d, 5910d], xsty=1
      cgoplot, lambda*1d10, spaxel(i,j,1,*)/I_lambda0, $
        background='white',color='red',axiscolor='black', charthick = 1, thick = 2
      cgoplot, lambda*1d10, spaxel(i,j,2,*)/I_lambda0,$
        background='white',color='blue',axiscolor='black', charthick = 1, thick = 2
      ; for transitions of two types of atoms with different veolocity distributions
      ; cgplot, v, P_v2, background='white', color='red', axiscolor='black'
       !p.multi=0
    ps_end
  endelse
  ;-------------------------------------------------------------------------------------
  ; save cache; debug
;  window, 2
;        cgplot, lambda*1d10, I1_lambda/I_lambda0, background='white',color='red',axiscolor='black', $
;        title = 'Relative intensity'
;      cgoplot, lambda*1d10, I2_lambda/I_lambda0, background='white',color='blue',axiscolor='black'
;      cgoplot, lambda*1d10, I0_lambda/I_lambda0, background='white',color='black',axiscolor='black'
;  save, /variables, filename='var.sav'
  return, [equivalent_width0, equivalent_width1, equivalent_width2, column_density(i,j)]
end

;    ; ------------ restrict range of lambda to avoid double peak due to negative velocities ----------
;    lambda_lowerlimit1 = min(v2lambda(v,lambda1))
;    lambda_upperlimit1 = max(v2lambda(v,lambda1))
;    print, 'limits of lambda as defined by the v_gaussian+v_bulk for Na I D1 is' $
;       + '[' + string(lambda_lowerlimit1) + ',' + string( lambda_upperlimit1) + ']'
;    lambda_lowerlimit2 = min(v2lambda(v,lambda2))
;    lambda_upperlimit2 = max(v2lambda(v,lambda2))
;    print, 'limits of lambda as defined by the v_gaussian+v_bulk for Na I D2 is' $
;       + '[' + string(lambda_lowerlimit2) + ',' + string( lambda_upperlimit2) + ']'
;    lambda_lowerlimit1_index = value_locate(lambda, lambda_lowerlimit1)
;    lambda_upperlimit1_index = value_locate(lambda, lambda_upperlimit1)
;    lambda_lowerlimit2_index = value_locate(lambda, lambda_lowerlimit2)
;    lambda_upperlimit2_index = value_locate(lambda, lambda_upperlimit2)
;    if lambda_lowerlimit1 lt lambda_min then begin
;      print, 'decrease lambda_min to display all data. execution halted'
;      stop      ; lambda_lowerlimit outside the lambda interval to be plotted on x-axis
;    endif
;    if lambda_upperlimit1 gt lambda_max then begin
;      print, 'increase lambda_max to display all data. execution halted'
;      stop      ; lambda_upperlimit on the right of the interval
;    endif
;    if lambda_lowerlimit2 lt lambda_min then begin
;      print, 'decrease lambda_min to display all data. execution halted'
;      stop      ; lambda_lowerlimit outside the lambda interval to be plotted on x-axis
;    endif
;    if lambda_upperlimit2 gt lambda_max then begin
;      print, 'increase lambda_max to display all data. execution halted'
;      stop      ; lambda_upperlimit on the right of the interval
;    endif
;    phi1(0:lambda_lowerlimit1_index)=0       ; null outside what's defined by the velocity distribution
;    phi1(lambda_upperlimit1_index: size(phi1, /n_elements)-1)=0
;    phi2(0:lambda_lowerlimit2_index)=0       ; null outside what's defined by the velocity distribution
;    phi2(lambda_upperlimit2_index: size(phi2, /n_elements)-1)=0