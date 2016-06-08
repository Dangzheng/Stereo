function depth = coarse2fine_dual(data, params,dual)
   level_cal = floor( ( log(params.minSize) - log(min(size(data.d0))) )/(log(params.scalefactor)) )+1;
   levels = min(params.level_max, round(level_cal));
   warps = dual.warps;
   for 
end