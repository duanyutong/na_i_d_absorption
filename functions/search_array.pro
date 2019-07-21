;+
; NAME:
; SEARCH_ARRAY
;
; PURPOSE:
;       Searches a one-dimensional array to find the first occurrence
;       of a value closest to a scalar value
;
; TYPE:
;       FUNCTION
;
; CATEGORY:
;       Arrays
;
; CALLING SEQUENCE:
; RESULT = SEARCH_ARRAY(INARR, TARGET)
;
; INPUTS:
;       INARR = the one-dimensional array to be searched.
; TARGET = the value to find in INARR (scalar or vector).
; 
; INPUT PARAMETERS:
;       NONE
; 
; KEYWORD PARAMETERS:
;       NONE
;
; OUTPUTS:
;       RESULT = index of INARR where INARR is closest to TARGET.
;        If TARGET is a vector, RESULT will be a vector of indices. 
;
; COMMON BLOCKS:
; NONE
;
; SIDE EFFECTS:
;       Finds only first occurance of TARGET.  If INARR is, for example, a
;       parabola and one desires to search beyond the maximum, pass in a
;       temporary array that is only the appropriate portion of INARR.
;
; RESTRICTIONS:
; NONE
;
; DEPENDENCIES:
; NONE
;
; MODIFICATION HISTORY:
; Written, Robert.Mallozzi@msfc.nasa.gov, October, 1994.
;-
; * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

FUNCTION SEARCH_ARRAY, INARR, TARGET

IF (N_ELEMENTS(TARGET) EQ 1) THEN BEGIN

   INDEX = WHERE(INARR EQ TARGET, CHECK)
   IF (CHECK EQ 0) THEN JUNK = MIN( ABS(INARR - TARGET), INDEX )

ENDIF ELSE BEGIN

   INDEX = LONARR(N_ELEMENTS(TARGET))
   FOR I=0, N_ELEMENTS(TARGET)-1 DO BEGIN
       JUNK = MIN( ABS(INARR - TARGET(I)), TEMP_INDEX )
       INDEX(I) = TEMP_INDEX
   ENDFOR

ENDELSE



RETURN, INDEX
END