pro cog
  ;+
  ; plot curve of growth
  ;-
  ; input
  density_coefficient_min = 1d-10
  density_coefficient_max = 1d10
  delta_x = 0.1
  @fc
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
  ;-----------------------------------------------------------------------------------
  x_min = alog10(density_coefficient_min)
  x_max = alog10(density_coefficient_max)
  x_index = dindgen(ceil((x_max - x_min)/delta_x))
  x = [x_min + x_index*delta_x, x_max]
  x_n = size(x, /n_elements)
  density_coefficient = exp(x)
  ;------------------------------------------------------------------------------------
  column_density = dblarr(x_n, /nozero)
  equivalent_width1 = dblarr(x_n, /nozero)
  equivalent_width2 = dblarr(x_n, /nozero)
  for i = 0, x_n-1, 1 do begin
    model1_return = dblarr(2, /nozero)
    model1_return = gas_model1_func(density_coefficient(i), /suppress_gfx)
    column_density(i) = model1_return(0)
    equivalent_width1(i) = model1_return(1)
    equivalent_width2(i) = model1_return(2)
  endfor
  
  ; plot curve of growth----------------------------------------------------------------
  ; w=0       ; window index
  window, w++
  cgplot, f1*column_density/lambda1, equivalent_width1/lambda1 $
    ,background='white',color='red',axiscolor='black', /xlog, /ylog
  window, w++
  cgplot, f2*column_density/lambda2, equivalent_width2/lambda2 $
    ,background='white',color='red',axiscolor='black', /xlog, /ylog
  stop
end