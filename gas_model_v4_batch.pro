windtype = 4
  for v_exp = 0,2 do begin
    gas_model_v4, windtype, v_exp, !pi/100*1
    gas_model_v4, windtype, v_exp, !pi/6*1
    gas_model_v4, windtype, v_exp, !pi/4*1
    gas_model_v4, windtype, v_exp, !pi/6*2
    gas_model_v4, windtype, v_exp, !pi/6*3
    gas_model_v4, windtype, v_exp, !pi/6*4
    gas_model_v4, windtype, v_exp, !pi/4*3
    gas_model_v4, windtype, v_exp, !pi/6*5
    gas_model_v4, windtype, v_exp, !pi/100*99
  endfor
end
