function b = baseline(G1, G2)

p1 = (-G1(1:3,1:3)')*(G1(1:3,4));
p2 = (-G2(1:3,1:3)')*(G2(1:3,4));
b = norm(p1 - p2,2);
end