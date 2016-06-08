function depth = l0norm(data, params)
   K = data.K;
   Iref = data.I1;
   I = data.I2;
   zeta = to_zeta(data.d);
   Grel = data.G;%左相机向右相机坐标转换的关系矩阵
   warp = warp_zeta(Iref, I, Grel, K, to_z(zeta));
   r = warp.r;
   a = warp.a;
   n = size(a,1)*size(a,2);
   
   zetak = zeta;
   grad = make_derivatives_2D_complete(size(z));
   Dx = grad.Kx;
   Dy = grad.Ky;
   
   
   Lnorm = make_linearOperator(size(zeta), K);
   lx = Lnorm.L_x;
   ly = Lnorm.L_y;
   lz = Lnorm.L_z;
   dxx = Dx*lx;
   dxy = Dx*ly;
   dxz = Dx*lz;
   dyx = Dy*lx;
   dyy = Dy*ly;
   dyz = Dy*lz;
   Dxx = (dxx)'*(dxx);
   Dxy = (dxy)'*(dxy);
   Dxz = (dxz)'*(dxz);
   Dyx = (dyx)'*(dxx);
   Dyy = (dyy)'*(dxy);
   Dyz = (dyz)'*(dxz);
   bot1 = Dxx + Dxy + Dxz + Dyx + Dyy + Dyz;
   while beta < betamax
       zeta_temp = zeta;
       %h-v subproblem
       zeta_flat = reshape(zeta.',[],1);
       zetak_flat = reshape(zetak.',[],1);
       x_norm = lx*zeta_flat;
       y_norm = ly*zeta_flat;
       z_norm = lz*zeta_flat;
       
       hpx = Dx*x_norm;
       hpy = Dx*y_norm;
       hpz = Dx*z_norm;
       
       vpx = Dy*x_norm;
       vpy = Dy*y_norm;
       vpz = Dy*z_norm;
       hv = hpx.^2 + hpy.^2 + hpz.^2 + vpx.^2 + vpy.^2 + vpz.^2;
       t = hv < lambd/beta;
       hpx(t) = 0;hpy(t) = 0;hpz(t) = 0;
       vpx(t) = 0;vpy(t) = 0;vpz(t) = 0;
       
       %zeta-subproblem
       
       a = spdiags(reshape(a.',[],1) ,0 ,n ,n);
       r = reshape(r.',[], 1);
       bot2 = a'*a;
       top1 = dxx*hpx + dxy*hpy + dxz*hpz + dyx*vpx +dyy*vpy + dyz*vpz;
       top2 = a'*a*zeta_flat - a'*zetak_flat;
       zeta = (bot1+bot2)\(top1+top2);
       zeta = reshape(full(zeta),fliplr(size(Iref)))';
   end
end
