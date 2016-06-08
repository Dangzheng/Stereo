function depth = l2norm(data, params)
   K = data.K;   
   Iref = data.I1;
   I = data.I2;
   z = data.d;
   Grel = data.G;%左相机向右相机坐标转换的关系矩阵
   warp = warp_z(Iref, I, Grel, K, z);
   r = warp.r;
   a = warp.a;
   n = size(a,1)*size(a,2);
   
   zk = z;
   grad = make_derivatives_2D_complete(size(z));
   Dx = grad.Kx;
   Dy = grad.Ky;
   Dx_x = Dx'*Dx;
   Dy_y = Dy'*Dy;
   
   itermax = params.itermax;
   beta = params.beta;
   energy = (reshape(r.',[],1)).^2 + beta*((Dx*(reshape(z.',[],1))).^2 + (Dy*(reshape(z.',[],1))).^2);
   %energy = (reshape(r.',[],1)+diag(reshape(a.',1,[]))*reshape((z-zk).',[],1)).^2 + beta*((Dx*(reshape(z.',[],1))).^2 + (Dy*(reshape(z.',[],1))).^2);
   energy = sum(sum(energy));
   %r_show = sum(sum(r.^2));
   energy_ = energy+1;
   disp(['---the one true energy: ',num2str(energy)]);
   for iter = 1:1:itermax
       z_temp = z;
       a = spdiags(reshape(a.',[],1),0,n,n);
       r = reshape(r.',[],1);
       zk = reshape(zk.',[],1);
       
       top = a'*a*zk - a'*r;
       top = sparse(top);
       bot = a'*a + beta*(Dx_x + Dy_y);
       z = bot\top;
       %z(z>params.maxz) = params.maxz;
       %z(z<params.minz) = params.minz;
       %energy1 = (r+a*(z-zk)).^2 + beta*((Dx*(reshape(z.',[],1))).^2 + (Dy*(reshape(z.',[],1))).^2);
       %energy1 = sum(sum(energy1));
       z = reshape(full(z),fliplr(size(Iref)))'; 
       
       %disp(['---after solve linear func: ',num2str(energy1)]);
       %imshow(uint8(z));
       zk = z_temp;
       warp = warp_z(Iref, I, Grel, K, z);
       r = warp.r;
       a = warp.a;
       %energy = (reshape(r.',[],1)+diag(reshape(a.',1,[]))*reshape((z-zk).',[],1)).^2 + beta*((Dx*(reshape(z.',[],1))).^2 + (Dy*(reshape(z.',[],1))).^2);
       energy = (reshape(r.',[],1)).^2 + beta*((Dx*(reshape(z.',[],1))).^2 + (Dy*(reshape(z.',[],1))).^2);
       energy = sum(sum(energy));
       
       if energy_ <= energy
           z = z_temp;
           disp(['the energy start to accelerate at ',num2str(iter),'th iter,the wrong energy is:',num2str(energy)]);
           break;
       end
       energy_ = energy;
       disp(['---energy (after update the parameters): ',num2str(energy)]);
   end
   %z(z>params.maxz) = params.maxz;
   %z(z<params.minz) = params.minz;
 
   hAxe = axes('Parent',gcf,... 
    'Units','pixels',...  
    'Position',[530 10 434 383]); 
   axes(hAxe);
   imshow(uint8(z*255/max(max(z))));
   pause(0.1)
   depth = z;

end
