 
 function depth = disp2depth(disp)
 disp = disp/8;
 [row,col] = size(disp);
 baseline = 2;
 Tx = - baseline;
 cx = col/2;
 cy = row/2;
 f = col;
 Q = [1,0,0,-cx;
     0,1,0,-cy;
     0,0,0,f;
     0,0,-1/Tx,0];
 [x,y] = size(disp);
 [x,y] = meshgrid(1:y,1:x);
 x = reshape(x.',1,[]);
 y = reshape(y.',1,[]);
 disp = reshape(disp.',1,[]);
 xyd1 = [x;y;disp;ones(1,row*col)];
 XYZW = Q*xyd1;
 xyz = XYZW(1:3,:)./repmat(XYZW(4,:),3,1);
 %figure(3);
 %plot3(xyz(1,:),xyz(2,:),-xyz(3,:),'.');
 %axis equal
 depth = reshape(xyz(3,:),fliplr([row,col]))';
 end