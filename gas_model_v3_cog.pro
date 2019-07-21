pro gas_model3_cog
  ;+
  ; plot curve of growth
  ;-
  ; input
  n_min = 1d-10
  n_max = 1d40
  delta_x = 2
  @fc
  ; --------------------------------- for transition 1-----------------------------------
  f1 = 0.640511d                         ; oscillator strength
  bkbj1 = 1d                      ; the ratio b_k/b_j is 1 for bound-bound transition
  lambda1 = 5889.950954d-10                    ; central wavelength for transition 1. 5890A in our case
  m_atomic1 = m_atomic_Na                  ; atomic mass of atom 1
  ; --------------------------------- for transition 2-----------------------------------
  f2 = 0.319913d                         ; oscillator strength
  bkbj2 = 1d                      ; the ratio b_k/b_j is 1 for bound-bound transition
  lambda2 = 5895.924237d-10                    ; central wavelength for transition 2. 5896A in our case
  m_atomic2 = m_atomic_Na                  ; atomic mass of atom 2
  ;-----------------------------------------------------------------------------------
  x_min = alog10(n_min)
  x_max = alog10(n_max)
  print, '['+string(x_min)+ ','+string(x_max)+ ']'
  x_index = dindgen(ceil((x_max - x_min)/delta_x))
  x = [x_min + x_index*delta_x, x_max]
  x_n = size(x, /n_elements)
  n = exp(x)
  ;------------------------------------------------------------------------------------
  column_density = dblarr(x_n, /nozero)
  equivalent_width0 = dblarr(x_n, /nozero)
  equivalent_width1 = dblarr(x_n, /nozero)
  equivalent_width2 = dblarr(x_n, /nozero)
  for i = 0, x_n-1, 1 do begin
    model_return = dblarr(2, /nozero)
    model_return = gas_model3_func(n(i), /suppress_gfx)
    equivalent_width0(i) = model_return(0)
    equivalent_width1(i) = model_return(1)
    equivalent_width2(i) = model_return(2)
    column_density(i) = model_return(3)
    print, 'finished'+string(i)+'of'+string(x_max)
  endfor
  
  ; plot curve of growth----------------------------------------------------------------
  ; w=0       ; window index
  ps_start, file='/users/duan/Documents/research/model3/cog.eps', /encap
    cgplot, alog10(f1*column_density/lambda1), alog10(equivalent_width1/lambda1) $
      ,background='white',color='red',axiscolor='black', charthick = 1, thick = 10 $
      ,title = 'Curve of Growth' $
      ,xtitle = 'log Nf/'+cgsymbol('lambda'), ytitle = 'log W/'+cgsymbol('lambda') $
      ,xr = [0, 21], xsty=1
    cgoplot, alog10(f2*column_density/lambda2), alog10(equivalent_width2/lambda2) $
      ,background='white',color='blue',axiscolor='black', charthick = 1, thick = 6
  ps_end
  stop
end