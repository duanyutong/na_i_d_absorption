function temp_spherical, r, theta
; temperature as a function of radius
; input: r, theta
common	shared_parameters, number_density_range, 					$
			number_density_coefficient, number_density_exp, 		$
			temp_coefficient, temp_exp, temp0,						$
			v_bulk_coefficient, v_bulk_exp, v_bulk0
			; share parameters across several functions
T = temp0 + temp_coefficient*r^temp_exp*(abs(cos(theta)))^0
return, T
end
