p1 = [0 0 ; 0 5; 8 0; 8 5; 1 4; 3 2; 5 4; 7 1]*21; % [mm]
p2 = [743 303; 536 430; 973 643; 777 784; 611 449; 752 484; 730 264; 906 628]; % [px]
n = size(p1,1);

A = zeros(2*n, 9);
for i = 1:n
    x = p1(i,1);
    y = p1(i,2);
    xp = p2(i,1);
    yp = p2(i,2);

    A(2*i-1,:) = [-x -y -1 0 0 0 x*xp y*xp xp];
    A(2*i,:) =  [0 0 0 -x -y -1 x*yp y*yp yp];
end

[U,S,V] = svd(A);
h = V(:,end);
H = reshape(h,[3 3])';

% Afficher l'image

img = imread('Points.jpeg'); 

figure; % Create a new figure window
imshow(img); % Display the image
hold on;

for i = 1:n
    xpt = H*[p1(i,:) 1]'*21;
    xpt = xpt/xpt(3);
    plot(xpt(1),xpt(2),'ro')
    hold on
end

