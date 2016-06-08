function Ir = imresize_m(img,sz) 

   if size(img) == sz
       Resize_image = img;
       return;
   end
   factors = sz./size(img);
   if any(factors < 1)%�ж�ͼ������С���ǷŴ���С�Ļ���Ҫ�Ӹ�˹ģ��
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
   %���������Ϊ��ֱ�Ӵ�python��ת�����������ԣ����������ǿ���������������Ȳ���һ��
   v = (v-1)*fy;v = v + (1/factors(1))/2 - 1 + 0.5;
   u_flat = reshape(u.',[],1);v_flat = reshape(v.',[],1);
   Ir = interp2(I_filter, u_flat,v_flat,'cubic');
   Ir = reshape(Ir,fliplr(sz))';
       
end