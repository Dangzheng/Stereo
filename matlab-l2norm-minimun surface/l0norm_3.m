function depth = l0norm_3(data, params)
   K = data.K;
   Iref = data.I1;
   I = data.I2;
   zeta = to_zeta(data.d);
   Grel = data.G;%左相机向右相机坐标转换的关系矩阵
   warp = warp_zeta(Iref, I, Grel, K, data.d);
   r = warp.r;
   a = warp.a;
   n = size(a,1)*size(a,2);
   
   zetak = zeta;

   Lnorm = make_linearOperator(size(zeta), K);
   lx = Lnorm.L_x;
   ly = Lnorm.L_y;
   lz = Lnorm.L_z;
   bot1 = lx'*lx + ly'*ly + lz'*lz;
   alpha = params.alpha;
   beta = params.beta;
   betamax = params.betamax;
   kappa = params.kappa;
   lambda = params.lambda;
   itermax = params.itermax;
   iter = 1;
   %while beta < betamax
   while iter < itermax
       zeta_tmp = zeta;
       %h-v subproblem
       zeta_flat = reshape(zeta.',[],1);
       zetak_flat = reshape(zetak.',[],1);
       hp = lx*zeta_flat;
       vp = ly*zeta_flat;
       dp = lz*zeta_flat;
       
       hvd = hp.^2 + vp.^2 + dp.^2;
       t = hvd < lambda/beta;
       hp(t) = 0;vp(t) = 0;dp(t) = 0;
       %zeta-subproblem
       a = spdiags(reshape(a.',[],1) ,0 ,n ,n);
       r = reshape(r.',[], 1);
       bot2 = a'*a;
       top1 = lx'*hp + ly'*vp + lz'*dp;
       top2 = a'*a*zetak_flat - a'*r;
       zetav = (beta*bot1 + alpha*bot2)\(beta*top1 + alpha*top2);
       zeta = reshape(zetav,fliplr(size(Iref)))';
       %z(z>params.maxz) = params.maxz;
       %zeta(zeta < params.minz) = to_zeta(params.minz);
       zeta(zeta < 0) = to_zeta(params.minz);
       warp = warp_zeta(Iref, I, Grel, K, to_z(zeta));
       r = warp.r;
       a = warp.a;
       energy = alpha*( reshape(r.',[],1) ).^2 +...
                beta*( (lx*zetav - hp).^2 + (ly*zetav - vp).^2 + (lz*zetav - dp).^2 );
       energy = sum(sum(energy)) + lambda*sum(sum(t));
       
       beta = beta * kappa;
       if iter == 1
           energy_ = energy + 1;
       end
       if energy_ <= energy
           zeta = zeta_tmp;
           disp(['the energy start to accelerate at ',num2str(iter),'th iter,the wrong energy is:',num2str(energy)]);
           break;
       end
       zetak = zeta_tmp;
       iter = iter + 1;
       disp(['---energy (after update the parameters): ',num2str(energy)]);
   end
   depth = to_z(zeta);
   
   hAxe = axes('Parent',gcf,... 
    'Units','pixels',...  
    'Position',[530 10 434 383]); 
   axes(hAxe);
   imshow(uint8(depth*255/max(max(depth))));
   pause(0.1)
end