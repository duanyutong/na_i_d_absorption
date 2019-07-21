function discrete_area, xdata, ydata

;+
; simple function to calculate area under curve
; defined by two discrete arrays
;-

  if size(xdata, /n_elements) ne size(ydata, /n_elements) then begin
    print, 'xdata and ydata dimensions mismatch. execution halted.'
    stop
  endif else begin
    n = size(xdata, /n_elements)
    xleft = xdata(0:n-2)
    xright = xdata(1:n-1)
    yleft = ydata(0:n-2)
    yright = ydata(1:n-1)
    b = xright - xleft
    h = (yleft + yright)/2
    area = total(b*h)
  return, area
  endelse
end
