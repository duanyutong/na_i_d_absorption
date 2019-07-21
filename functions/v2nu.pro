function v2nu, v, nu0

  @fc                   ; load funcamental constants

	beta = v/c
	nu=nu0*sqrt((1-beta)/(1+beta))
	return, nu
	
end
