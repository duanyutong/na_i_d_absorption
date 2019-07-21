function v_add, v0, u
  ;+
  ; accurate velocity addition according special relativity
  ; input parameters
  ; v0: original velocity in rest frame
  ; u: velocity of moving frame with respect to rest frame
  ;-
  
  @fc                   ; load funcamental constants

  v = (u + v0)/(1d + u*v0/c^2d)
  return, v
  
end
