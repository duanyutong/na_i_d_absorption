function v_bulk_spherical, r, theta
; bulk velocity of gas as a function of radius
; input: r, theta
; v - r^2
common	shared_parameters, number_density_range, 					$
			number_density_coefficient, number_density_exp, 		$
			temp_coefficient, temp_exp, temp0,						$
			v_bulk_coefficient, v_bulk_exp, v_bulk0
			; share parameters across several functions
v_bulk = v_bulk0 + v_bulk_coefficient*r^(v_bulk_exp)
return, v_bulk
end
