; docformat = 'rst'
;+
; This is an example program to demontrate how to create an image plot with
; contours overlaid on it with Coyote Graphics routines.
;
; :Categories:
;    Graphics
;    
; :Examples:
;    Save the program as "image_plot_with_contours.pro" and run it like this::
;       IDL> .RUN image_plot_with_contours
;       
; :Author:
;    FANNING SOFTWARE CONSULTING::
;       David W. Fanning 
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: david@idlcoyote.com
;       Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 23 January 2013 by David W. Fanning.
;
; :Copyright:
;     Copyright (c) 2013, Fanning Software Consulting, Inc.
;-
PRO Image_Plot_with_Contours
   restore, file='/users/duan/Documents/research/model3/variables.sav' ; restore session 
   ; Example Gaussian data.
   image = equivalent_width(*,*,0)
   
   ; Set up variables for the contour plot. Normally, these values 
   ; would be passed into the program as positional and keyword parameters.
   minValue = Min(image)
   maxValue = Max(image)
   nLevels = 10
   xtitle = 'X Axis'
   ytitle = 'Y Axis'
   position =   [0.125, 0.125, 0.9, 0.800]
   cbposition = [0.125, 0.865, 0.9, 0.895]
   cbTitle = 'Data Value'
   
   ; Set up a "window" for the plot. The PostScript output will have
   ; the same aspect ratio as the graphics window on the display.
   cgDisplay, 600, 600, Title='Image Plot with Contours'
   
   ; Set up colors for contour plot.
   cgLoadCT, 33
   
   ; Display the image on the display. Keep its aspect ratio.
   cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, $
       /Axes, XTitle=xtitle, YTitle=ytitle, Position=position, $
       XRange=[-grid_radius, grid_radius], YRange=[-grid_radius, grid_radius], /Keep_Aspect
   
   ; Overplot the contours.
   contourLevels = cgConLevels(image, NLevels=10, MinValue=minValue)
   cgContour, image, Levels=contourLevels, /OnImage, Color='charcoal'
   
   ; Draw the color bar. Fit it to the location of the image.
   cgColorbar, Position=cbposition, Range=[MinValue, MaxValue], $
      Title=cbTitle, /Fit

END ;*****************************************************************

; This main program shows how to call the program and produce
; various types of output.

  ; Display the plot in a graphics window.
  Image_Plot_with_Contours
  
  ; Display the plot in a resizeable graphics window.
  cgWindow, 'Image_Plot_with_Contours', WXSize=600, WYSize=650, $
      WTitle='Image Plot with Contours in Resizeable Graphics Window'
  
  ; Create a PostScript file.
  PS_Start, 'image_plot_with_contours.ps'
  Image_Plot_with_Contours
  PS_End
  
  ; Create a PNG file with a width of 600 pixels.
  cgPS2Raster, 'image_plot_with_contours.ps', /Portrait, /PNG, Width=600

END