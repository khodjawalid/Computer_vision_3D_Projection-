function [H] =  myDLT(m1,m2)
%H = [0 0 0 ,-m2',m1(2)*m2'; m2' ,0 0 0,-m1(1)*m2';-m1(2)*m2',m1(1)*m2', 0 0 0]

% m2 : points world
% m1 : points images

format long
taille = size(m1)
nbPoint = taille(1)

a = zeros(3,3)
A = zeros(2*nbPoint, 9)

for i =1:nbPoint
M1 = [m1(i,:)]'   % point image
M2 = [m2(i,:)]'   % point world

a = [0 0 0 ,-M2',M1(2)*M2'; M2' ,0 0 0,-M1(1)*M2';-M1(2)*M2',M1(1)*M2', 0 0 0]

if i == 1 
    A(i,:) = a(1,:)
    A(i+1,:) = a(2,:)
else
    A(2*i-1,:) = a(1,:)
    A(2*i,:) = a(2,:)
end

[U,S,V] = svd(A)


H = V(:,end);

end