  ps_start, file='/users/duan/documents/research/model3/colormap^' $
                    +strtrim(number_density_exp,2)+'_v^' $
                    +strtrim(v_bulk_exp,2)+'_t^' $
                    +strtrim(temp_exp,2)+'.eps', /encap
    !p.multi=[0, 3, 1, 0, 0]
    margin = [3, 3, 3, 3]
   ; --------------------- 1st graph equivalent width -------------------
   ; Set up variables for the contour plot. Normally, these values 
   ; would be passed into the program as positional and keyword parameters.
   IMAGE = equivalent_width(*,*,0)*1D9
   MINVALUE = MIN(IMAGE)
   MAXVALUE = MAX(IMAGE)
   nLevels = 10
   xtitle = 'X Axis'
   ytitle = 'Y Axis'
   position =   [0.125, 0.125, 0.9, 0.800]
   cbposition = [0.125, 0.865, 0.9, 0.895]
   cbTitle = 'Data Value/nm'
   
   ; Set up a "window" for the plot. The PostScript output will have
   ; the same aspect ratio as the graphics window on the display.
   ; cgDisplay, 600, 600, Title='Image Plot with Contours'
   
   ; Set up colors for contour plot.
   cgLoadCT, 33
   
   ; Display the image on the display. Keep its aspect ratio.
   cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, $
       /Axes, XTitle=xtitle, YTitle=ytitle, Position=position, $
       XRange=[-grid_radius, grid_radius], YRange=[-grid_radius, grid_radius], /Keep_Aspect, title='Equivalent Width'
   
   ; Overplot the contours.
   contourLevels = cgConLevels(image, NLevels=10, MinValue=minValue)
   cgContour, image, Levels=contourLevels, /OnImage, Color='charcoal'
   
   ; Draw the color bar. Fit it to the location of the image.
   cgColorbar, Position=cbposition, Range=[MinValue, MaxValue], $
      Title=cbTitle, /Fit
   ; ----------------------- 2nd graph column density ----------------------------
   ; Set up variables for the contour plot. Normally, these values 
   ; would be passed into the program as positional and keyword parameters.
   image = column_density/1d14
   minValue = Min(image)
   maxValue = Max(image)
   nLevels = 10
   xtitle = 'X Axis'
   ytitle = 'Y Axis'
   position =   [0.125, 0.125, 0.9, 0.800]
   cbposition = [0.125, 0.865, 0.9, 0.895]
   cbTitle = 'Data Value/1E14'
   
   ; Set up a "window" for the plot. The PostScript output will have
   ; the same aspect ratio as the graphics window on the display.
   ; cgDisplay, 600, 600, Title='Image Plot with Contours'
   
   ; Set up colors for contour plot.
   cgLoadCT, 33
   
   ; Display the image on the display. Keep its aspect ratio.
   cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, $
       /Axes, XTitle=xtitle, YTitle=ytitle, Position=position, $
       XRange=[-grid_radius, grid_radius], YRange=[-grid_radius, grid_radius], /Keep_Aspect, title='Column Density'
   
   ; Overplot the contours.
   contourLevels = cgConLevels(image, NLevels=10, MinValue=minValue)
   cgContour, image, Levels=contourLevels, /OnImage, Color='charcoal'
   
   ; Draw the color bar. Fit it to the location of the image.
   cgColorbar, Position=cbposition, Range=[MinValue, MaxValue], $
      Title=cbTitle, /Fit
   ; ------------------------------- 3rd graph galactic disc and wind ----------------------
   ; Set up variables for the contour plot. Normally, these values 
   ; would be passed into the program as positional and keyword parameters.
   image = transpose(grid.plane(grid_radius,*,*))
   minValue = Min(image)
   maxValue = Max(image)
   nLevels = 10
   xtitle = 'Z Axis (Line of Sight)'
   ytitle = 'Y Axis'
   position =   [0.125, 0.125, 0.9, 0.800]
   cbposition = [0.125, 0.865, 0.9, 0.895]
   cbTitle = 'Data Value'
   
   ; Set up a "window" for the plot. The PostScript output will have
   ; the same aspect ratio as the graphics window on the display.
   ; cgDisplay, 600, 600, Title='Image Plot with Contours'
   
   ; Set up colors for contour plot.
   cgLoadCT, 33
   
   ; Display the image on the display. Keep its aspect ratio.
   cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, $
       /Axes, XTitle=xtitle, YTitle=ytitle, Position=position, $
       XRange=[-grid_radius, grid_radius], YRange=[-grid_radius, grid_radius],$
       /Keep_Aspect, title='Plane Orientation'
   PARTVELVEC, grid.v_bulk_z(grid_radius,*,*), grid.v_bulk_y(grid_radius,*,*), $
                grid.z(grid_radius,*,*), grid.y(grid_radius,*,*), $
                /over, fraction=1
;   ; Overplot the contours.
;   contourLevels = cgConLevels(image, NLevels=10, MinValue=minValue)
;   cgContour, image, Levels=contourLevels, /OnImage, Color='charcoal'
;   
   ; Draw the color bar. Fit it to the location of the image.
   cgColorbar, Position=cbposition, Range=[MinValue, MaxValue], $
      Title=cbTitle, /Fit
  !p.multi=0
  ps_end