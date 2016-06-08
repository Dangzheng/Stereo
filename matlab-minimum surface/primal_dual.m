function depth = primal_dual(data, params)
   K = data.K;
   Iref = data.I1;
   I = data.I2;
   Grel = data.G;
   Lnorm = make_linearOperator(size(Iref), K);
   L_normal = Lnorm.L;
   d = data.d;
   p = data.p;
   for k = 1:1:params.warps
       disp(['---warp: ',num2str(k),' ---maxZ = ',num2str(params.maxz)]);
       warp = warp_zeta(Iref, I, Grel, K, data.d);
       solve = solve_area_pd(warp, d, p, params, L_normal);
       d = solve.d;
       p = solve.p;
   end 
   depth = solve;
end