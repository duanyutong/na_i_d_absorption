pro gas_model_v3_colormaps
  ;+
  ; plot colormaps
  ;-
  ; INPUT
  ;     spaxel(x, y, trasition, spectrum)
  ;     equivalent_width(x, y, transition)
  ;     column_density(x,y)
  @fc
  ;------------------------------ main -------------------------------
  restore, file='/users/duan/Documents/research/model3/variables.sav' ; restore session
  ps_start, file='/users/duan/documents/research/model3/colormap_wt^' $
                      +strtrim(floor(windtype),2)+'_res' $
                      +strtrim(floor(grid_radius)*2,2)+'_incl' $
                      +strtrim(inclination,2) $
                      , /encap
    !p.multi=[0, 3, 1, 0, 0]
    margin = [3, 3, 3, 3]
   ; --------------------- 1st graph equivalent width -------------------
   ; Set up variables for the contour plot. Normally, these values 
   ; would be passed into the program as positional and keyword parameters.
   image = equivalent_width(*,*,0)*1D9
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
       /Axes, XTitle='X Distance / kpc', YTitle='Y Distance / kpc', Position=position, $
       XRange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       YRange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       /Keep_Aspect, title='Equivalent Width / nm'
   ; Draw the color bar. Fit it to the location of the image.
   cgColorbar, Position=cbposition, Range=[min(image), max(image)], /Fit
   ; ----------------------- 2nd graph column density ----------------------------
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
       /axes, xtitle='X Distance / kpc', ytitle='Y Distance / kpc', position=position, $
       xrange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       yrange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       /keep_aspect, title='Column Density / 1E14 '+textoidl('m^{-2}')
   ; draw the color bar. fit it to the location of the image.
   cgcolorbar, position=cbposition, range=[min(image), max(image)], /fit
   ; ------------------------------- 3rd graph galactic disc and wind ----------------------
   ; Set up variables for the contour plot. Normally, these values 
   ; would be passed into the program as positional and keyword parameters.
   velocity_magnitude = sqrt(grid.v_bulk_z(grid_radius,*,*)^2 + grid.v_bulk_y(grid_radius,*,*)^2)
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
   cgimage, image+plane, stretch=1, minvalue=min(image), maxvalue=max(image), $
       /axes, xtitle='Z Distance (Line of Sight) / kpc', ytitle='Y Distance / kpc', position=position, $
       xrange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       yrange=[-r_galaxy/pc/1d3, r_galaxy/pc/1d3], $
       /keep_aspect, title='Velocity Field / '+textoidl('km\cdot')+textoidl('s^{-1}')
   PARTVELVEC, grid.v_bulk_z(grid_radius,*,*), grid.v_bulk_y(grid_radius,*,*), $
                r_galaxy/pc/1d3/grid_radius*grid.z(grid_radius,*,*), $
                r_galaxy/pc/1d3/grid_radius*grid.y(grid_radius,*,*), $
                /over, fraction=0.2
   ; Draw the color bar. Fit it to the location of the image.
   cgColorbar, Position=cbposition, Range=[min(image), max(image)], /Fit
   !p.multi=0
  ps_end
  stop
end
