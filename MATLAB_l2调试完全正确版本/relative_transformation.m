function Grel = relative_transformation(Gref,g)
%Ϊʲô�ܾ���������������������ϵת��
Grel = zeros(3,4);
Grel(1:3,1:3) = g(1:3,1:3)*Gref(1:3,1:3)';
Grel(1:3,4) = g(1:3,4) - Grel(1:3,1:3)*Gref(1:3,4);
end
