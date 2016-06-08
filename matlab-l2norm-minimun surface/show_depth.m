
function show_depth(depth,i)
 [x,y] = size(depth);
 [x,y] = meshgrid(1:y,1:x);
 x = reshape(x.',1,[]);
 y = reshape(y.',1,[]);
 depth = reshape(depth.',1,[]);
 xyz(1,:) = x;xyz(2,:) = y;
 xyz(3,:) = depth;
 figure(i);
 plot3(xyz(1,:),xyz(2,:),xyz(3,:),'.');
 axis equal;
end