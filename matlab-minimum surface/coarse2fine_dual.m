function depth = coarse2fine_dual(data, dual)
   level_cal = floor( ( log(dual.minSize) - log(min(size(data.d0))) )/(log(dual.scalefactor)) )+1;
   levels = min(dual.levels, round(level_cal));
   %计算primal-dual所需要的参数
   warps = dual.warps;
   iteraions = dual.iterations;
   
   L = dual.lambda;
   lambdual = zeros(levels,1);
   for i = 1:1:size(lambdual,1);
       lambdual(i) = L*( (1/dual.scalefactor).^i );
   end
   dim_dual = 3;
   
   for i = levels:-1:dual.stop_level
      %跟据新的图片的分辨率计算相应的内参矩阵
      level_sz = round(size(data.d0)*((dual.scalefactor)^(i-1)));
      factor = level_sz./size(data.d0);
      factor(1,3) = 1;
      K = data.K .*((repmat(factor,3,1))');
      disp(['---level ',num2str(i),' ---size:(', num2str(level_sz(1)),',',num2str(level_sz(2)),')']);
      if i == levels
          d = imresize(data.d0,level_sz,'bicubic');
          p = zeros([size(d), dim_dual]);
      else
          d = imresize(d, level_sz,'bicubic');
          ptmp = p;
          p = zeros([size(d), dim_dual]);
          for j = 1:1:dim_dual
              p(:,:,j) = imresize(ptmp(:,:,j),level_sz, 'bicubic');
          end
      end
      data_scaled.I1 = imresize(data.I1,level_sz,'bicubic');
      data_scaled.I2 = imresize(data.I2,level_sz,'bicubic');
      data_scaled.G = data.G;
      data_scaled.K = K;
      data_scaled.d = d;
      data_scaled.p = p;
      dual.lambda = lambdual(i);
      depth = primal_dual(data_scaled, dual);
      data_scaled.d = depth.d;
      p = depth.p;
      hAxe = axes('Parent',gcf,... 
      'Units','pixels',...  
      'Position',[530 10 434 383]); 
      axes(hAxe);
      imshow(uint8(depth.d*255/max(max(depth.d))));
      pause(0.0001);
   end
end