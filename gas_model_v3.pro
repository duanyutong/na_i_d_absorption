; for best view set tab size to 4 characters. in vim run :set tabstop=4
pro gas_model_v3, windtype_input, inclination_input, SUPPRESS_GFX=suppress_gfx
	@fc		; load funcamental constants
	common	shared_parameters, number_density_range, 				$
			number_density_coefficient, number_density_exp, 		$
			temp_coefficient, temp_exp, temp0,						$
			v_bulk_coefficient, v_bulk_exp, v_bulk0
			; share parameters across several functions
			; DO NOT CHANGE as order of these variables MATTERS as IDL is dumb
; =================== define input variables =====================
; variables for MAXWELLIAN distribution : m_atomic, T
; variables for guassian distribution   : mu, sigma
; Geometry assumptions:	observing along z direction
;			major axis of the ellipse observed always
;			is the x-axis
;			galactic center is at (0,0) on x-y plane
; --------------------- model options ---------------------
; spot = [0, 0]		; location of the pixel observed on x-y plane
windtype = windtype_input	; 1, sperical
							; 2, filled biconical
							; 3, unfilled biconical
; windtype = 2d		
; grid_radius = 10d			; number of pixels along radius
;							; half length of each dimension
vdt = 2						; velocity distribution type:	1 for maxwellian
;											 				2 for gaussian
distance_offset = 1d	; distance from galactic plane. for plotting purpose
bicone_angle = !pi/6d	; angle between cone and galactic plane (theta = pi/2)
bicone_thickness = !pi/13d	; for windtype 3
; --------------------- observables ------------------------
; r_proj = 0.5*r_galaxy				; projected radius. measured vertically
inclination = inclination_input		; we can measure original & projected
									; area or radius
; inclination = !pi*3/4				; inclination angle between galactic
									; disc and z+ axis. range [0, pi)
; used !pi/2 (vertical) and !pi*3/4 (45 degrees) for senior presentation
r_galaxy = 10d3*pc		; actual radius of the galaxy. measured horizontally
r_wind = r_galaxy*0.7	; radius of the gas
center_perturb = r_galaxy/grid_radius/2	; amount of perturbation for center
										; location to calculate quantities
v_hubble = 200d3						; velocity due to hubble expansion
										; typically 200km/s (M82)
; --------------------- dependences ----------------------
sigma = 100d3	; standard deviation for gaussian velocity distribution
mu = 0d			; mean for gaussian
number_density_range = [3d-7, 1d-8]		; from center to outer boundary
number_density_exp =  -2				; [-3, 0]
temp_exp =  -1							; [-3, 0]
temp_range = [1d4,1d3]					; for maxwell velocity distribution
v_bulk_exp = 1							; [ 0, 3]
v_bulk_range = [50d3,500d3]
; ----------------- plotting settings ------------------- 
; without considering velocity due to hubble expansion
lambda_min = 5800d-10		; lower limit of line profile domain in the 
							; plot and all calculations
lambda_max = 5950d-10		; upper limit of line profile domain
delta_lambda = 0.10d-10		; increment of wavelength
v_gaussian_max = 800d3		; upper limit of gaussian velocity
delta_v_gaussian = 0.05d3	; increment of thermal velocity
I0 = 1d						; original intensity of light beam
I_lambda0 = I0				; use a function insted for a radiation 
							; profile that depends on wavelength.
							; eg. blackbody 
; ------------ Na I D1 transition parameters -----------
f1 = 0.319913d				; oscillator strength
bkbj1 = 1d					; b_k/b_j is 1 for bound-bound transition
lambda1 = 5897.56661715d-10	; central wavelength for transition 1,
							; 5890A in our case
m_atomic1 = m_atomic_Na		; atomic mass of atom 1
; -------------- Na I D2 ------------------------
f2 = 0.640511d				; oscillator strength
bkbj2 = 1d					; b_k/b_j is 1 for bound-bound transition
lambda2 = 5891.58326415d-10	; central wavelength for transition 2,
							; 5896A in our case
m_atomic2 = m_atomic_Na		; atomic mass of atom 2
; ============= set coefficients of dependence relations  =================
; exponential form of functions assumed. given boundary values, calculate
; coefficients in functions
number_density_coefficient =	1/r_galaxy*								$
								(										$
									number_density_range(1)				$
									/									$
									number_density_range(0)				$
								)										$
								^(1/number_density_exp)
temp_coefficient =	(temp_range(1)-temp_range(0))						$
					/													$
					(r_galaxy^temp_exp - center_perturb^temp_exp)
v_bulk_coefficient = 	(v_bulk_range(1)-v_bulk_range(0))				$
						/												$
						(r_galaxy^v_bulk_exp - center_perturb^temp_exp)
temp0 = temp_range(1) - temp_coefficient * r_galaxy^temp_exp
v_bulk0 = v_bulk_range(1) - v_bulk_coefficient * r_galaxy^v_bulk_exp
; ====================== construct lattice =====================
; ---------------------- basic stuff  ------------------
side = 2*grid_radius+1		; total number of pixels along each dimension
xhat = [1, 0, 0]			; unit vectors in cartesian space
yhat = [0, 1, 0]
zhat = [0, 0, 1]
normal =	[							$	; normal vector
				0,						$	; with respect to galactic disc
				sin(inclination-!pi/2),	$	; always in the y-z plane
				cos(inclination-!pi/2)	$
			]
plane_co1 = normal(0)		; coefficients to specify galactic plane
plane_co2 = normal(1)
plane_co3 = normal(2)
; ----------------------- allocate memory ------------------
occupy = dblarr(side,side,side)				; just an empty template
grid =	{								$
			x:				occupy,		$	; index of x coordinate
			y:				occupy,		$	; index of y coordinate
			z:				occupy, 	$	; index of z coordinate
			r:				occupy,		$	; radial coordinate
			theta:			occupy, 	$	; polar angle with respect to 
										$	; the normal direction
										$	; which is the polar axis
			number_density:	occupy,		$
			temp:			occupy, 	$
			v_bulk:			occupy,		$
			v_bulk_x:		occupy,		$
			v_bulk_y:		occupy,		$
			v_bulk_z:		occupy,		$
			plane:			occupy  	$	; as logical matrix. value is 1
										$	; if the plane lies in the block
		}
; --------------------- loop for each cell in lattice ----------------
; assign index number with respect to central origin, (i-grid_radius), to
; every slice
for i = 0, side-1 do begin
	occupy = replicate(i-grid_radius,side,side)
	grid.x(i,*,*) = occupy
	grid.y(*,i,*) = occupy
	grid.z(*,*,i) = occupy		
endfor
for i = 0, side-1 do begin
	for j = 0, side-1 do begin
		for k = 0, side-1 do begin	
			position =	[						$
							grid.x(i,j,k),		$
							grid.y(i,j,k),		$
							grid.z(i,j,k)		$
						]			; dimensionless index numbers
									; in the cartesian grid
        	grid.r(i,j,k) = norm(position)*r_galaxy/grid_radius
									; radial distance
        	grid.theta(i,j,k) = acos	(                   $ 
                                			total(position*normal)$
                            			    /(norm(position)*norm(normal)) 
										)	; angle between two vectors
			position_spherical = [grid.r(i,j,k), grid.theta(i,j,k)]
			; assume azimuthally symmetric wind
        	; ----------------- outside the wind ----------------
	        if grid.r(i,j,k) gt r_wind then begin
				grid.number_density(i,j,k) = 0d		; set density to 0
													; and we are done
				; ~~~~~~~ determine if block is a plane block ~~~~~~~ 
				; los_logic = grid.theta(i,j,*) gt !pi/2
				; plane_index = SEARCH_ARRAY(los_logic, 0)
				; plane_index = plane_index(0)
				; grid.plane(i,j,plane_index) = 1
				distance = abs	(								$
									plane_co1*grid.x(i,j,k)+	$
									plane_co2*grid.y(i,j,k)+	$
									plane_co3*grid.z(i,j,k)		$
								)
                          / sqrt(plane_co1^2+plane_co2^2+plane_co3^2)
						  ; distance to the galactic disk
				if distance le distance_offset then begin
					grid.plane(i,j,k)=1
				endif	$
				else begin
					grid.plane(i,j,k)=0
				endelse
			endif	$
        	; ----------------- inside the wind --------------------
			else begin
				; ~~~~~~~ determine if block is a plane block ~~~~~~~ 
				; los_logic = grid.theta(i,j,*) gt !pi/2
				; plane_index = SEARCH_ARRAY(los_logic, 0)
				; plane_index = plane_index(0)
				; grid.plane(i,j,plane_index) = 1
				distance = abs	(								$
									plane_co1*grid.x(i,j,k)+	$
									plane_co2*grid.y(i,j,k)+	$
									plane_co3*grid.z(i,j,k)		$
								)
                          / sqrt(plane_co1^2+plane_co2^2+plane_co3^2)
						  ; distance to the galactic disk
				if distance le distance_offset then begin
					grid.plane(i,j,k)=1
				endif	$
				else begin
					grid.plane(i,j,k)=0
				endelse
				; ~~~~~~~~ fill in physical parameters ~~~~~~~~
				case windtype of
					1:	begin			; spherial wind
                  			grid.number_density(i,j,k) = $
								number_density_spherical(position_spherical)
							grid.temp(i,j,k) = $
								temp_spherical(position_spherical)
							v_bulk =	position / norm(position) *			$
										abs	(								$
												v_bulk_spherical	(		$
														position_spherical	$
																	)		$
											)
										; the function v_bulk_spherical 
										; returns velocities from given
										; spherical coordinates (r, theta)
										; in cartesian vector form
							grid.v_bulk(i,j,k) = norm(v_bulk)
							grid.v_bulk_x(i,j,k) = v_bulk(0)
							grid.v_bulk_y(i,j,k) = v_bulk(1)
							grid.v_bulk_z(i,j,k) = v_bulk(2)
						end
					2:	begin			; filled biconical wind
							if abs(grid.theta(i,j,k)-!pi/2) lt bicone_angle	$
								then begin
								grid.number_density(i,j,k)=0d
							endif	$
							else begin
								grid.number_density(i,j,k) =				$
									number_density_spherical	(			$
														position_spherical	$
																)
								grid.temp(i,j,k) = temp_spherical	(		$
														position_spherical	$
																	)
								v_bulk =	position / norm(position) *		$
											abs	(							$
													v_bulk_spherical	(	$
														position_spherical	$
																		)	$
												)
											; the function v_bulk_spherical 
											; returns velocities from given
											; spherical coordinates (r, theta)
											; in cartesian vector form
								grid.v_bulk(i,j,k) = norm(v_bulk)
								grid.v_bulk_x(i,j,k) = v_bulk(0)
								grid.v_bulk_y(i,j,k) = v_bulk(1)
								grid.v_bulk_z(i,j,k) = v_bulk(2)
							endelse
						end
					3:	begin				; unfilled biconical wind
							if												$
								(abs(grid.theta(i,j,k)-bicone_angle) le 	$
									bicone_thickness						$
								)											$
								||											$
								(abs(grid.theta(i,j,k)-(!pi-bicone_angle)	$
									 le bicone_thickness					$
								)											$
							then begin		; on either side of galactic disk
								grid.number_density(i,j,k) =				$
									number_density_spherical(				$
														position_spherical	$
															)				$
								grid.temp(i,j,k) =							$
									temp_spherical(position_spherical)
								v_bulk = position/norm(position)*			$
									abs(v_bulk_spherical(position_spherical))
									; the function v_bulk_spherical 
									; returns velocities from given
									; spherical coordinates (r, theta)
									; in cartesian vector form
									; velocities are relative to galaxy
								grid.v_bulk(i,j,k) = norm(v_bulk)
								grid.v_bulk_x(i,j,k) = v_bulk(0)
								grid.v_bulk_y(i,j,k) = v_bulk(1)
								grid.v_bulk_z(i,j,k) = v_bulk(2)
							endif	$
							else begin
								grid.number_density(i,j,k)=0d
							endelse
						end			; finish windtype 3
				endcase				; finish windtype cases
			endelse					; finish inside the wind condition
		endfor
	endfor
    print, 'finished contructing plane x = ' + string(i)
endfor
; -------------- fill galactic center by perturbation ----------
; as some functions yield null values when r=0 or theta=0
position = [0,center_perturb,0]    ;  index for shifted core position
grid.r(grid_radius,grid_radius,grid_radius) = center_perturb
grid.theta(grid_radius,grid_radius,grid_radius) =	$
	acos	(										$
				total(position*normal)				$	; dot product
				/									$
				(norm(position)*norm(normal))		$	; product of norms
			)
grid.plane(grid_radius,grid_radius,grid_radius) = 1
position_spherical = [center_perturb, 1d-3]
grid.number_density(grid_radius,grid_radius,grid_radius) = 			$
	number_density_spherical(position_spherical)
grid.temp(grid_radius,grid_radius,grid_radius) =					$
	temp_spherical(position_spherical)
v_bulk = position/norm(position)*abs(v_bulk_spherical(position_spherical))    
grid.v_bulk(grid_radius,grid_radius,grid_radius) = norm(v_bulk)
grid.v_bulk_x(grid_radius,grid_radius,grid_radius) = v_bulk(0)
grid.v_bulk_y(grid_radius,grid_radius,grid_radius) = v_bulk(1)
grid.v_bulk_z(grid_radius,grid_radius,grid_radius) = v_bulk(2)
; ========================= calculate spaxels ========================
; ---------------- create empty matrices --------------
; create wavelength array
lambda_index = dindgen(ceil((lambda_max - lambda_min)/delta_lambda))
lambda = [lambda_min + lambda_index*delta_lambda, lambda_max]
lambda_n = size(lambda, /n_elements)
; first two dimensions for pixel position, third dimension for transition
; selection, forth dimension for storing spectrum
spaxel = dblarr(side,side,3,lambda_n)		; spaxel
tau = dblarr(side,side,3,lambda_n)			; optical depth
column_density = dblarr(side,side)
equivalent_width = dblarr(side,side,3)
v_gaussian_index = dindgen(ceil((2*v_gaussian_max/delta_v_gaussian)))
v_gaussian = 	[															$
					-v_gaussian_max + v_gaussian_index*delta_v_gaussian,	$
					v_gaussian_max											$
				]							; gaussian spread
; ------------ integrate along line of sight ----------- 
for i = 0, side-1 do begin
	for j = 0, side-1 do begin			; loop over all spaxels
		; i = spot(0)+grid_radius		; for a specific spaxel
		; j = spot(1)+grid_radius
		; if line of sight passes through galactic disk plane=1 for some cells
		; the sum of plane is greater than 0
		if total(grid.plane(i,j,*)) gt 0 then begin
			for k = side-1, 0, -1 do begin	; integrate along -z direction
				; stop when hitting galactic plane
				if grid.plane(i,j,k) eq 1 then break
					; ~~~~~~ compute integrated absoption coefficient ~~~~~~~~
					s1 = acs_int(f1, bkbj1, lambda1, grid.temp(i,j,k))
					s2 = acs_int(f2, bkbj2, lambda2, grid.temp(i,j,k))
					; ~~~~~~~~ determine wavelengths of interest ~~~~~~~~~~
					; info_lambda1_min = min(v2lambda(v,lambda1))
					; info_lambda2_min = min(v2lambda(v,lambda2))
					; info_lambda1_max = max(v2lambda(v,lambda1))
					; info_lambda2_max = max(v2lambda(v,lambda2))
					; ~~~~~~ velocity and wavelength distribution ~~~~~~~~~
					; array of absolute velocities relative to observer
					v = v_add	(									$
									grid.v_bulk_z(i,j,k),			$
									v_add(v_gaussian, v_hubble)		$
								)
					; v_n = size(v, /n_elements)
					; target wavelength space to velocity space
					; convert to the grid cell rest frame
					; for the first transition
					v_rest1 = v_add	(
										lambda2v(lambda, lambda1),			$
										v_add	(							$
													grid.v_bulk_z(i,j,k),	$
													-v_hubble				$
												)							$
									)
					; for the second transition
					v_rest2 = v_add	(
										lambda2v(lambda, lambda2),			$
										v_add	(							$
													grid.v_bulk_z(i,j,k),	$
													-v_hubble				$
												)							$
									)
					case vdt of
						1:	begin		; maxwell distribution
								T = grid.temp(i,j,k)
								phi1 =										$
									P_v2phi									$
										(									$
											v2lambda(v_rest1, lambda1),		$
											lambda1,						$
											ATOMICMASS=m_atomic1,			$
											TEMP=T							$
										)
								phi2 =										$
									P_v2phi									$
										(									$
											v2lambda(v_rest2, lambda2),		$
											lambda2,						$
											ATOMICMASS=m_atomic2,			$
											TEMP=T							$
										)
							end
						2:  begin      ; gaussian distribution
								phi1 = P_v2phi	(							$
													v2lambda	(			$
														v_rest1, lambda1	$
																),			$
													lambda1,				$
													MEAN=mu,				$
													SD=sigma				$
												)	$
                                )   ; assume same gaussian for all locations
								phi2 = P_v2phi	(							$
													v2lambda	(			$
														v_rest2, lambda2	$
																),			$
													lambda1,				$
													MEAN=mu,				$
													SD=sigma				$
												)	$
                end
          endcase
          ; ~~~~ absorption cross section for each specific wavelength ~~~~~~
          
          s1_lambda = s1*phi1
          s2_lambda = s2*phi2
          s0_lambda = s1_lambda + s2_lambda
          ; column density
          column_density_block = grid.number_density(i,j,k)*$
			  r_galaxy/grid_radius
          ; optical depth
          tau(i,j,0,*) = tau(i,j,0,*) + s0_lambda*column_density_block
          tau(i,j,1,*) = tau(i,j,1,*) + s1_lambda*column_density_block
          tau(i,j,2,*) = tau(i,j,2,*) + s2_lambda*column_density_block
          column_density(i,j) = column_density(i,j) + column_density_block
        endfor
      endif
      ; ---------------------- resulting intensity ---------------------- 
      spaxel(i,j,0,*) = I_lambda0*exp(-tau(i,j,0,*))
      spaxel(i,j,1,*) = I_lambda0*exp(-tau(i,j,1,*))
      spaxel(i,j,2,*) = I_lambda0*exp(-tau(i,j,2,*))
      ; ------------------------ equivalent width -----------------------
      equivalent_width(i,j,0) = discrete_area$
		  (lambda, 1-spaxel(i,j,0,*)/I_lambda0)
      equivalent_width(i,j,1) = discrete_area$
		  (lambda, 1-spaxel(i,j,1,*)/I_lambda0)
      equivalent_width(i,j,2) = discrete_area$
		  (lambda, 1-spaxel(i,j,2,*)/I_lambda0)
    endfor
    print, 'finished integrating over spaxal column x = ' + string(i)
  endfor
  ; save cache; debug
  save, /all, /comm, /variables, $
	  file='/users/duan/Documents/research/model3/variables_wt'+$
      strtrim(floor(windtype),2)+'_res'+$
      strtrim(floor(grid_radius)*2,2)+'_incl'+$
      strtrim(inclination,2)+$
      '.sav'
                                            
  ; ========================== plots =========================
  
  if keyword_set(suppress_gfx) then begin
  endif else begin
  
    ; ------------------------- spaxel_display ------------------------- 
    
    ps_start, file='/users/duan/documents/research/model3/spaxel_wt'$
                        +strtrim(floor(windtype),2)+'_res' $
                        +strtrim(floor(grid_radius)*2,2)+'_incl' $
                        +strtrim(inclination,2) $
                        +'.eps', /encap, /nomatch
      x = grid_radius+1     ; set center coordinates
      y = grid_radius+1
      spec0 = spaxel(x,y,0,*)
      spec1 = spaxel(x,y,1,*)
      spec2 = spaxel(x,y,2,*)
      margin = [2, 2, 2, 2]
      cgplot,   lambda*1d10, spec0/I_lambda0, $
		  background='white',color='black',axiscolor='black', $
        xtitle = 'Wavelength / '+cgsymbol('Angstrom'), $
		ytitle = 'Relative Flux' $
        , xr = [5870d, 5930d], xsty=1, charthick = 1, thick = 2
      cgoplot,  lambda*1d10, spec1/I_lambda0, $
		  background='white',color='red',axiscolor='black', $
		  charthick = 1, thick = 2
      cgoplot,  lambda*1d10, spec2/I_lambda0, $
		  background='white',color='blue',axiscolor='black', $
		  charthick = 1, thick = 2
    ps_end
    ; colormaps
    ps_start, file='/users/duan/documents/research/model3/colormap_wt^' $
                      +strtrim(floor(windtype),2)+'_res' $
                      +strtrim(floor(grid_radius)*2,2)+'_incl' $
                      +strtrim(inclination,2) $
                      +'.eps', /encap
    !p.multi=[0, 2, 3, 0, 0]
    margin = [2, 2, 1, 0]
    ;
    ; |----------------------------------|
    ; | n-r       W_eq        v_mean  |
    ; | v-r       CD              P.O.          |
    ; |----------------------------------|
    ; --------------------- 1st graph v-r
    ; 
    ; 
    ; 
    ; --------------------- 1st graph equivalent width -------------------
    ; 
    ; Set up variables for the contour plot. Normally, these values 
    ; would be passed into the program as positional and keyword parameters.
    image = equivalent_width(*,*,0)*1D9  ; convert m to nm
    nLevels = 10
    position =   [0.125, 0.125, 0.9, 0.800]
    cbposition = [0.125, 0.865, 0.9, 0.895]
   
   ; Set up a "window" for the plot. The PostScript output will have
   ; the same aspect ratio as the graphics window on the display.
   ; cgDisplay, 600, 600, Title='Image Plot with Contours'
   
   ; Set up colors for contour plot.
   cgLoadCT, 33
   
   ; Display the image on the display. Keep its aspect ratio.
   cgImage, image, Stretch=1, MinValue=MIN(IMAGE), MaxValue=MAX(IMAGE), $
       /Axes, XTitle='X Distance / kpc', YTitle='Y Distance / kpc', $
	   Position=position, $
       XRange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       YRange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       /Keep_Aspect, title='Equivalent Width / nm'
   ; Draw the color bar. Fit it to the location of the image.
   cgColorbar, Position=cbposition, Range=[min(image), max(image)], /Fit
   
   ; ----------------------- 2nd graph column density ---------------------
   ; 
   ; Set up variables for the contour plot. Normally, these values 
   ; would be passed into the program as positional and keyword parameters.
   image = column_density/1d14
   nLevels = 10
   position =   [0.125, 0.125, 0.9, 0.800]
   cbposition = [0.125, 0.865, 0.9, 0.895]
   
   ; Set up a "window" for the plot. The PostScript output will have
   ; the same aspect ratio as the graphics window on the display.
   ; cgDisplay, 600, 600, Title='Image Plot with Contours'
   
   ; Set up colors for contour plot.
   cgLoadCT, 33
   
   ; Display the image on the display. Keep its aspect ratio.
   cgimage, image, stretch=1, minvalue=min(image), maxvalue=max(image), $
       /axes, xtitle='X Distance / kpc', ytitle='Y Distance / kpc', $
	   position=position, $
       xrange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       yrange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       /keep_aspect, title='Column Density / 1E14 '+textoidl('m^{-2}')
   ; draw the color bar. fit it to the location of the image.
   cgcolorbar, position=cbposition, range=[min(image), max(image)], /fit
   
   ; ------------------------------- 3rd graph galactic disc and wind -------
   ; 
   ; Set up variables for the contour plot. Normally, these values 
   ; would be passed into the program as positional and keyword parameters.
   velocity_magnitude = sqrt$
	   (grid.v_bulk_z(grid_radius,*,*)^2 + grid.v_bulk_y(grid_radius,*,*)^2)
   image = transpose(velocity_magnitude/1d3)
   plane = transpose(grid.plane(grid_radius,*,*))*max(image)
   nLevels = 10
   position =   [0.125, 0.125, 0.9, 0.800]
   cbposition = [0.125, 0.865, 0.9, 0.895]
   
   ; Set up a "window" for the plot. The PostScript output will have
   ; the same aspect ratio as the graphics window on the display.
   ; cgDisplay, 600, 600, Title='Image Plot with Contours'
   
   ; Set up colors for contour plot.
   cgLoadCT, 33
   
   ; Display the image on the display. Keep its aspect ratio.
   cgimage, image+plane, stretch=1,$
	   minvalue=min(image), maxvalue=max(image), $
       /axes, xtitle='Z Distance (Line of Sight) / kpc', $
	   ytitle='Y Distance / kpc', position=position, $
       xrange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       yrange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       /keep_aspect, $
	   title='Velocity Field / '+textoidl('km\cdot')+textoidl('s^{-1}')
   partvelvec, grid.v_bulk_z(grid_radius,*,*), $
	   grid.v_bulk_y(grid_radius,*,*), $
                r_galaxy/pc/1d3/grid_radius*grid.z(grid_radius,*,*), $
                r_galaxy/pc/1d3/grid_radius*grid.y(grid_radius,*,*), $
                /over, fraction=0.3
   ; Draw the color bar. Fit it to the location of the image.
   cgColorbar, Position=cbposition, Range=[min(image), max(image)], /Fit
   !p.multi=0
  ps_end
;  r_index = dindgen(ceil((r_max - r_min)/delta_r))
;  r_array = [r_min + r_index*delta_r, r_max]
;  r_n = size(r, /n_elements)
;    w=0       ; window index
;    window, w++     ; velocity distribution of Na
;    i = spot(0)+grid_radius
;    j = spot(1)+grid_radius
;    z = grid.z(i,j,*)*r_galaxy/grid_radius
  endelse

stop
end


