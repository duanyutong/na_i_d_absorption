function acs_int, f, bkbj, nu, T

;calculate integrated absorption cross section
;for a given upward oscillator strength
;
;bkbj     : the ratio b_k/b_j is 1 for bound-bound transition
;f        : upward oscillator strength
;lambda   : central wavelength for a certain transition
;T        : temperature

@fc									; load funcamental constants
s_u = (!pi*e^2)/(m_e*c)*f/epsilon	; uncorrected for stimulated emission
s = s_u*(1-bkbj*exp(-(h*nu)/(k*T)))	; corrected for stimulated emission
return, s

end
