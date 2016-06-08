function Grel = relative_transformation(Gref,g)
%为什么总觉得是右相机向左相机坐标系转？
Grel = zeros(3,4);
Grel(1:3,1:3) = g(1:3,1:3)*Gref(1:3,1:3)';
Grel(1:3,4) = g(1:3,4) - Grel(1:3,1:3)*Gref(1:3,4);
end
