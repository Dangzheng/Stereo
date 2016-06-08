function warp = warp_zeta(Iref, I, Grel, K, z)
[x,y] = size(z);
   [x,y] = meshgrid(1:y,1:x);
   fx = K(1,1);fy = K(2,2);
   cx = K(1,3);cy = K(2,3);
   Xn = ones(3, numel(x));
   x_flat = reshape(x.',1,[]);%x,y,z 是按照行进行展开的
   y_flat = reshape(y.',1,[]);
   z_flat = reshape(z.',1,[]);
   Xn(1,:) = (x_flat - cx)/fx;%因为python是从0开始计数，所以这里有些不一样
   Xn(2,:) = (y_flat - cy)/fy;
   
   X = Xn;
   X(4,:) = 1./z_flat;
   [row,col] = size(Grel);
   if row < 4
       Grel(4,:) = [0,0,0,1];
   end
   X2 = Grel * X;%X2为右相机坐标系下的坐标表示
   x2 = K * X2(1:3,:);
   x2 = x2./repmat(x2(3,:),3,1); %x2也是安照行进行排列的
   I = double(I);
   Iw = interp2(I,x2(1,:),x2(2,:),'cubic');
   Iw = reshape(Iw,fliplr(size(I)))';
   Iw(isnan(Iw)) = 0;
   warp.r = Iw - Iref;%residential
   warp.r(isnan(warp.r)) = 0;
   %imshow(uint8(Iw))
   
   
   dX = Grel(1:3,1:3)*Xn./z_flat;%derivation of X2 wrt. z
   x_pos = x2(1:2,:);
   x_pos(1,:) = x_pos(1,:) + 0.5*ones(1,length(x_pos));
   Ix_pos = interp2(I,x_pos(1,:),x_pos(2,:),'cubic');
   x_neg = x2(1:2,:);
   x_neg(1,:) = x_neg(1,:) - 0.5*ones(1,length(x_neg));
   Ix_neg = interp2(I,x_neg(1,:),x_neg(2,:),'cubic');
   
   y_pos = x2(1:2,:);
   y_pos(2,:) = y_pos(2,:) + 0.5*ones(1,length(y_pos));
   Iy_pos = interp2(I,y_pos(1,:),y_pos(2,:),'cubic');
   y_neg = x2(1:2,:);
   y_neg(2,:) = y_neg(2,:) - 0.5*ones(1,length(y_neg));
   Iy_neg = interp2(I,y_neg(1,:),y_neg(2,:),'cubic');
   
   gIwx = Ix_pos - Ix_neg;
   gIwy = Iy_pos - Iy_neg;

   X2 = X2./repmat(X2(4,:),4,1); %dehomogenize
   
   z2 = X2(3,:).^2;
   dT = zeros(2,length(X2));
   dT(1,:) = fx*ones(1,length(dX))./X2(3,:).*dX(1,:) - fx*ones(1,length(dX)).*X2(1,:)./z2.*dX(3,:);
   dT(2,:) = fy*ones(1,length(dX))./X2(3,:).*dX(2,:) - fy*ones(1,length(dX)).*X2(2,:)./z2.*dX(3,:);
   %一直到dT都没有进行过转换，所以一直都是按招行进行排列的
   gIwx_flat = reshape(gIwx.',1,[]);
   gIwy_flat = reshape(gIwy.',1,[]);
   warp.a = reshape(gIwx_flat.*dT(1,:) + gIwy_flat.*dT(2,:),fliplr(size(z)))';
   warp.a(isnan(warp.a)) = 0;
end