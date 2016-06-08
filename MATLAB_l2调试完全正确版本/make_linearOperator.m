function Lnorm = make_linearOperator(size_zeta, K)
   m = size_zeta(1);n = size_zeta(2);
   fx = K(1,1);fy = K(2,2);
   [x,y] = meshgrid(1:n,1:m);
   x_flat = reshape(x.',1,[]);%x,y,z 是按照行进行展开的
   y_flat = reshape(y.',1,[]);
   x_hat = (x_flat - cx)/fx;%因为python是从0开始计数，所以这里有些不一样
   y_hat = (y_flat - cy)/fy;
   
   grad = make_derivatives_2D_complete(size_zeta);
   Kx = grad.Kx;Ky = grad.Ky;
   
   spId = speye(m*n);
   spXhat = spdiags(reshape(x_hat.',[],1),0,numel(x_hat),numel(x_hat));
   spYhat = spdiags(reshape(y_hat.',[],1),0,numel(x_hat),numel(x_hat));
   
   Lnorm.L_x = -Kx/fy;
   Lnorm.L_y = -Ky/fx;
   Lnorm.L_z = spXhat*Kx/fy +spYhat*Ky/fx + 2*spId/(fx*fy);
   
end