function solve = solve_area_pd(warp, d, p, params, L)
   check = params.check;
   epsilon = params.epsilon;
   lambda = params.lambda;
   It = warp.r;
   Ig = warp.a;
   %Iw = warp.Iw;
   
   zeta0 = to_zeta(d);
   zeta = zeta0;
   zeta_ = zeta0;
   
   %tau = zeros(size(L,2));%compute preconditioners;
   %sigma = zeros(size(L,1));
   tau = reshape(1./sum(abs(L),1),[],1);
   sigma = reshape(1./sum(abs(L),2),[],1);
   tau(isinf(tau)) = 0;
   sigma(isinf(sigma)) = 0;
   
   for k = 1:1:params.iterations
      p = flat(p);
      p = p + sigma.* ( L * reshape(zeta_.',[],1) );%dual update
      p = reshape(p,[],3)';
      
      normp = sqrt(p(1,:).^2 + p(2,:).^2 + p(3,:).^2);
      normp(normp <= 1) = 1;%reprojection 
      p = p./repmat(normp,3,1);
      p = reshape(p.',[],1);
      
      zeta_ = zeta;%remember zeta
      zeta = reshape(zeta.',[],1);
      
      zeta = zeta - tau.*(L'*p);
      zeta = reshape(zeta, fliplr(size(zeta0)))';
      tau = reshape(tau, fliplr(size(zeta0)))';
      r = It + Ig.*(zeta - zeta0);
      th = epsilon*ones(size(tau)) + lambda.*tau.*(Ig.^2);
      idx1 = r > th;
      idx2 = r < -th;
      idx3 = abs(r) <= th;
      
      zeta(idx1) = zeta(idx1) - lambda.*tau(idx1).*Ig(idx1);
      zeta(idx2) = zeta(idx2) + lambda.*tau(idx2).*Ig(idx2);
      zeta(idx3) = (zeta(idx3) - lambda.*tau(idx3).*Ig(idx3).* ...
      (It(idx3) - zeta0(idx3).*Ig(idx3))/epsilon)./ ...
      (1 + lambda.*tau(idx3).*Ig(idx3).^2/epsilon);
      tau = reshape(tau.',[], 1);
      
      zeta(zeta < params.minz) = params.minz;
      zeta(zeta > params.maxz) = params.maxz;
      zeta_ = 2*zeta - zeta_;
      grad = L*reshape(zeta.',[],1);
      grad = reshape(grad,3,[]);
      normgrad = sqrt(grad(1,:).^2 + grad(2,:).^2 + grad(3,:).^2);
      
      r = huber_function(It + (zeta - zeta0).* Ig, epsilon);
      energy = sum(sum( reshape(normgrad,fliplr(size(r)))' + lambda*r));
      if mod(k,check) == 0
          disp(['iter: ',num2str(k),' energy: ',num2str(energy)]);
      end
   p = reshape(p,[],3);
   p1(:,:,1) = reshape(p(:,1),fliplr(size(zeta0)))';
   p1(:,:,2) = reshape(p(:,2),fliplr(size(zeta0)))';
   p1(:,:,3) = reshape(p(:,3),fliplr(size(zeta0)))';
   p = p1;
   end
   solve.d = to_z(zeta);
   solve.p = p;
   
end