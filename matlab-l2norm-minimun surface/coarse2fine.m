function depth = coarse2fine(data, params)
%计算所需要的level的数目
level_cal = floor( ( log(params.minSize) - log(min(size(data.d0))) )/(log(params.scalefactor)) )+1;
levels = min(params.level_max, round(level_cal));

for i = levels:-1:params.stop_level
   %跟据新的图片的分辨率计算相应的内参矩阵
   level_sz = round(size(data.d0)*((params.scalefactor)^(i-1)));
   factor = level_sz./size(data.d0);
   factor(1,3) = 1;
   K = data.K .*((repmat(factor,3,1))');
   disp(['---level ',num2str(i),' ---size:(', num2str(level_sz(1)),',',num2str(level_sz(2)),')']);
   if i == levels
       d = imresize(data.d0,level_sz,'bicubic');
   else
       d = imresize(d, level_sz,'bicubic');
   end
   data_scaled.I1 = imresize(data.I1,level_sz,'bicubic');
   %imshow(uint8(data_scaled.I1));
   data_scaled.I2 = imresize(data.I2,level_sz,'bicubic');
   data_scaled.G = data.G;
   data_scaled.K = K;
   data_scaled.d = d;
   d = l2norm(data_scaled,params);
   %d = l2_6norm(data_scaled, params);
end
depth = d;
end