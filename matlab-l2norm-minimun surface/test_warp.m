   
   for z = 1:512
   x = 1,y = 1
   c = 434,r = 383;
   K = [c 0 c/2;0 c r/2;0 0 1];
   fx = K(1,1);fy = K(2,2);
   cx = K(1,3);cy = K(2,3);
   Xn = ones(3, numel(x));
   Grel = [1 0 0 -baseline;0 1 0 0;0 0 1 0];
   Xn(1,:) = (x - cx)/fx;%��Ϊpython�Ǵ�0��ʼ����������������Щ��һ��
   Xn(2,:) = (y - cy)/fy;
   
   X = Xn;
   X(4,:) = 1/z;
   [row,~] = size(Grel);
   if row < 4
       Grel(4,:) = [0,0,0,1];
   end
   X2 = Grel * X;%X2Ϊ���������ϵ�µ������ʾ
   x2 = K * X2(1:3,:);
   x2 = x2./repmat(x2(3,:),3,1); %x2Ҳ�ǰ����н������е�
   test(:,z) = x2;
   end