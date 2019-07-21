function number_density_spherical, r, theta
; number of absorbing atoms in unit volume along the radial direction
; input: r, theta
; rho - r^-2
common	shared_parameters, number_density_range, 					$
			number_density_coefficient, number_density_exp, 		$
			temp_coefficient, temp_exp, temp0,						$
			v_bulk_coefficient, v_bulk_exp, v_bulk0
			; share parameters across several functions
number_density =	number_density_coefficient*						$
					r^(number_density_exp)							$
					*(abs(cos(theta)))^0
return, number_density
end
