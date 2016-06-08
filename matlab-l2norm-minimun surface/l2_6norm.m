function depth = l2_6norm(data, params)
   K = data.K;
   Iref = data.I1;
   I = data.I2;
   zeta = to_zeta(data.d);
   Grel = data.G;%左相机向右相机坐标转换的关系矩阵
   warp = warp_zeta(Iref, I, Grel, K, data.d);
   r = huber_function(warp.r,params.epsilon);
   a = warp.a;
   n = size(a,1)*size(a,2);
   
   zetak = zeta;
   grad = make_derivatives_2D_complete(size(zeta));
   Dx = grad.Kx;
   Dy = grad.Ky;
   
   Lnorm = make_linearOperator(size(zeta), K);
   lx = Lnorm.L_x;
   ly = Lnorm.L_y;
   lz = Lnorm.L_z;
   
   %预计算L与梯度算子的乘积
   dxx = Dx*lx;
   dxy = Dx*ly;
   dxz = Dx*lz;
   dyx = Dy*lx;
   dyy = Dy*ly;
   dyz = Dy*lz;
   %预计算用于计算导数为零时要使用的梯度算子
   Dxx = (dxx)'*(dxx);
   Dxy = (dxy)'*(dxy);
   Dxz = (dxz)'*(dxz);
   Dyx = (dyx)'*(dxx);
   Dyy = (dyy)'*(dxy);
   Dyz = (dyz)'*(dxz);
   
   bot1 = Dxx + Dxy + Dxz + Dyx + Dyy + Dyz;
   lambda = params.lambda;
   iter = 1;
   itermax = params.itermax;
   while iter < itermax
      zeta_tmp = zeta;
      a = spdiags(reshape(a.',[],1),0,n,n);
      r = reshape(r.',[],1);
      zeta_flat = reshape(zeta.',[],1);
      zetak_flat = reshape(zetak.',[],1);
      top = -a'*r + a'*a*zetak_flat;
      bot2 = a'*a;
      zetav = (lambda*bot1 + bot2)\top;
      zetav(zetav > to_zeta(params.maxz)) = to_zeta(params.maxz);
      zetav(zetav < 0) = to_zeta(params.minz);
      zeta = reshape(full(zetav),fliplr(size(Iref)))';
      warp = warp_zeta(Iref, I, Grel, K, to_z(zeta));
      r = huber_function(warp.r,params.epsilon);
      a = warp.a;
      energy = ( reshape(r.',[],1) ).^2 +...
          (dxx*zetav).^2 + (dxy*zetav).^2 + (dxz*zetav).^2 + ...
          (dyx*zetav).^2 + (dyy*zetav).^2 + (dyz*zetav).^2 ;
      energy = sum(sum(energy));
      if iter == 1
         energy_ = energy + 1;
      end
       
      if energy_ <= energy
         zeta = zeta_tmp;
         disp(['the energy start to accelerate at ',num2str(iter),'th iter,the wrong energy is:',num2str(energy)]);
         break;
      end      
      iter = iter + 1;
      zetak = zeta_tmp;
   end
    depth = to_z(zeta);
    hAxe = axes('Parent',gcf,... 
    'Units','pixels',...  
    'Position',[530 10 434 383]); 
    axes(hAxe);
    imshow(uint8(depth*255/max(max(depth))));
    pause(0.01)
end
