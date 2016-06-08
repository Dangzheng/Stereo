function Ir = imresize_m(img,sz) 

   if size(img) == sz
       Resize_image = img;
       return;
   end
   factors = sz./size(img);
   if any(factors < 1)%判断图像是缩小还是放大，缩小的话就要加高斯模糊
       sigmas = (1./factors)/3.0;
       sigmas = mean(sigmas);
       gaussian = fspecial('gaussian',[3 3],sigmas);
       I_filter =  imfilter(img,gaussian,'replicate');
   else
       I_filter = img;
   end
   [u,v] = meshgrid(1:sz(2),1:sz(1));
   [row, col] = size(img);
   fx = col/sz(2);fy = row/sz(1);
   u = (u-1)*fx;u = u + (1/factors(2))/2 - 1 + 0.5;% sample from correct position
   %这个步骤因为是直接从python上转换过来的所以，坐标索引是可能有问题的所以先测试一下
   v = (v-1)*fy;v = v + (1/factors(1))/2 - 1 + 0.5;
   u_flat = reshape(u.',[],1);v_flat = reshape(v.',[],1);
   Ir = interp2(I_filter, u_flat,v_flat,'cubic');
   Ir = reshape(Ir,fliplr(sz))';
       
end