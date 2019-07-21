pro spaxel_display, x, y
  @fc                                 ; load funcamental constants
  common shared_parameters, number_density_coefficient, number_density_exp, number_density0, $
                            temp_coefficient, temp_exp, temp0, $
                            v_bulk_coefficient, v_bulk_exp, v_bulk0          ; share parameters with functions
  
  restore, file='/users/duan/Documents/research/model3/variables.sav' ; restore session
  ps_start, file='/users/duan/documents/research/model3/spaxel_wt'$
                      +strtrim(floor(windtype),2)+'_res' $
                      +strtrim(floor(grid_radius)*2,2)+'_incl' $
                      +strtrim(inclination,2) $
                      +'.eps', /encap
    spec0 = spaxel(x,y,0,*)
    spec1 = spaxel(x,y,1,*)
    spec2 = spaxel(x,y,2,*)
    window
    cgplot,   lambda*1d10, spec0/I_lambda0, background='white',color='black',axiscolor='black', $
      xtitle = 'Wavelength / '+cgsymbol('Angstrom'), ytitle = 'Relative Flux' $
      , xr = [5870d, 5930d], xsty=1, charthick = 1, thick = 2
    cgoplot,  lambda*1d10, spec1/I_lambda0, background='white',color='red',axiscolor='black'
    cgoplot,  lambda*1d10, spec2/I_lambda0, background='white',color='blue',axiscolor='black'
  ps_end
  stop
end