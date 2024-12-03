%chargement de l'image dans la matrice 


% for n=1:23
%   images{n} = imread(sprintf('Right\\right_frame%02d.png',n));
% end
% 
% imshow(images{22})
    
h = zeros(3,3)
m1= [1 2 1; 2 3 1; 1 1 1; 2 4 1]
m2= [1 5 1; 2 1 1; 5 4 1; 8 4 1]

H= myDLT(m1,m2)

h(1,1:3) = H(1:3)
h(2,1:3) = H(4:6)
h(3,1:3) = H(7:9)
format short

m1_h = h*[1 5 1]'
m1_h_n = [m1_h(1)/m1_h(3) m1_h(2)/m1_h(3) 1]

m2_h = inv(h)* [1 2 1]'
m2_h_n = [m2_h(1)/m2_h(3) m2_h(2)/m2_h(3) 1]

% Img0= imread('Right/right_frame00.png');
% figure()
% imshow(Img0);
% 
% Img0= imread('Right/right_frame00.png');
% figure(2)
% imshow(Img0);

% K = cameraParams.Intrinsics.IntrinsicMatrix'
% t1 = cameraParams.TranslationVectors(1,:)'
% R1 = cameraParams.RotationMatrices(:,:,1)'
% P1 = K*[R1,t1]
% 
% point12 = cameraParams.WorldPoints(12,:)
% point = [point12 0 1]'
% 
% 
% ppt_caclcule = P1*point
% ppt_calcule_normalise = [ppt_caclcule(1)/ppt_caclcule(3) ; ppt_caclcule(2)/ppt_caclcule(3)]
% 
% pro_point12 = cameraParams.ReprojectedPoints(12,:,1)'

% K =
% 
%    1.0e+03 *
% 
%     1.3620         0    0.7391
%          0    1.3609    0.9940
%          0         0    0.0010




