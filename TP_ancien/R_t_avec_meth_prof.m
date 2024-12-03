% 
% K = cameraParams.Intrinsics.IntrinsicMatrix'
% 
% % point World
% point2 = [cameraParams.WorldPoints(2,:) 1]
% point24 = [cameraParams.WorldPoints(24,:) 1]
% point38 = [cameraParams.WorldPoints(38,:) 1]
% point10 = [cameraParams.WorldPoints(10,:) 1]
% 
% % point images
% pro_point2 = cameraParams.ReprojectedPoints(2,:,1)
% pro_point24 = cameraParams.ReprojectedPoints(24,:,1)
% pro_point38 = cameraParams.ReprojectedPoints(38,:,1)
% pro_point10 = cameraParams.ReprojectedPoints(10,:,1)
% 
% % ensemble de point world
% M0 = [point2;point24;point38;point10]
% % ensemble de point image
% M1 = [pro_point2;pro_point24;pro_point38;pro_point10]
% 
% % calcul de l'homographie 9*1
% 
% % faut mettre M1 les point images et M0 les point world
% H1= myDLT(M1,M0)
% 
% % matrice d'homographie 3*3
% h1(1,1:3) = H1(1:3);
% h1(2,1:3) = H1(4:6);
% h1(3,1:3) = H1(7:9);
% 
% k_1H = inv(K)*h1
% alpha_r1 = k_1H(:,1)
% alpha_r2 = k_1H(:,2)
% alpha_t  = k_1H(:,3)
% 
% alpha2_r3 = cross(alpha_r1,alpha_r2)
% 
% alpha4 = det([alpha_r1 alpha_r2 alpha2_r3])
% 
% alpha_ = -(alpha4)^(1/4)
% alpha = (alpha4)^(1/4)
% % calcul des vecteurs de la matrice de rotation
% 
% r1_ = alpha_r1 /alpha_
% r2_ = alpha_r2/alpha_
% r3_ = alpha2_r3 / (alpha_)^2
% 
% 
% r1 = alpha_r1 /alpha
% r2 = alpha_r2/alpha
% r3 = alpha2_r3 / (alpha)^2
% 
% %format long
% 
% 
% % matrice de rotation calculée
% % R orienté suivant -z
% R_ = [r1_ r2_ r3_]
% % R orienté suivant z
% R = [r1 r2 r3]
% 
% % matrice de rotation donnée par matlab
% Rm =  cameraParams.RotationMatrices(:,:,1)'
% % vecteur de translation calculé
% 
% % vecteur de translation
% % suivant -z
% t_ = alpha_t /alpha_
% % suivant z
% t = alpha_t /alpha
% 
% % vecteur de translation donné par matlab
% tm = cameraParams.TranslationVectors(1,:)'
% 
% 
% P = K*[R_ t_]
% 
% Pm = K*[Rm tm]
% 
% % test 
% format short
% N = P*[cameraParams.WorldPoints(2,:) 0 1; cameraParams.WorldPoints(24,:) 0 1;cameraParams.WorldPoints(38,:) 0 1;cameraParams.WorldPoints(10,:) 0 1]'
% 
% p1 = N(:,1)
% p1n = p1(1:2)/p1(3)
% pro_point2'
% 
% 
% p2 = N(:,2)
% p2n = p2(1:2)/p2(3)
% pro_point24'
% 
% p3 = N(:,3)
% p3n = p3(1:2)/p3(3)
% pro_point38'
% 
% p4 = N(:,4)
% p4n = p4(1:2)/p4(3)
% pro_point10'



% % récuperer les images de la video
%  F = extractFrames("projet_fond_noir.mp4")



% avec fond noir :
% detection carré bleu :
% figure
% imagesc(I(:,:,1)<40 &  I(:,:,2)>0 & I(:,:,2) < 130 & I(:,:,3) > 120)
% % detection carré rouge :
% %figure
% imagesc(I(:,:,1)>180 & I(:,:,2)>25 & I(:,:,2) < 150 &  I(:,:,3) > 55  & I(:,:,3) < 170)
% % detection du carré vert :
% %figure
% imagesc(I(:,:,1)< 15 &  I(:,:,2)>140 & I(:,:,3)>120 & I(:,:,3)<220 )
% % detection du carré Jaune :
% %figure
% imagesc(I(:,:,1)> 190 &  I(:,:,2)> 190  & I(:,:,3)>100 & I(:,:,3)<210)



vectb = [];
vectr = [];
vectg = [];
vecty = [];

% for n=001:435
%   %images{n}
%   I = imread(sprintf('myVideoFrames2/000%03d.jpg',n));
%   
%   
% % caclul du centre de gravité du carré bleu
% Db = I(:,:,1)<40 &  I(:,:,2)>0 & I(:,:,2) < 130 & I(:,:,3) > 120;
% Db=double(Db);
% [Lb,numb]=bwlabel(Db,4); %etiquetage région 4 connexité , L contient les regions  
% %numbegionb=numb; % affichage nombre de région
% Hb=hist(Lb(:),[0:numb]); 
%     for j=2:length(Hb)
%         if Hb(j)<110 % on a choisi le seuil 20 pour determiner les petites region
%             Lb(Lb==j-1)=0; %suppression des petite regions -- L==j-1 récuprère l etiquete j-1
%         else
%             [yb,xb]=find(Lb==j-1);
%             vectb=[vectb; [mean(xb),mean(yb)]]; %calcul du centre de gravité des région
%                 % c est une concaténation vect=[vect; [mean(x),mean(y)]];
%         end 
%     end
%     
% %     figure(8),hold on, imagesc(Lb),colormap,title('image finale avec le centre de gravité');
% %     plot(vectb(:,1),vectb(:,2), '+b'); % trajectoire centre de gravité avec des carré bleu
%    
% % caclul du centre de gravité du carré rouge
%     
% Dr = I(:,:,1)>180 & I(:,:,2)>25 & I(:,:,2) < 150 &  I(:,:,3) > 55  & I(:,:,3) < 170;
% Dr=double(Dr);
% [Lr,numr]=bwlabel(Dr,4); %etiquetage région 4 connexité , L contient les regions  
% %numregionr=numr; % affichage nombre de région
% Hr=hist(Lr(:),[0:numr]); 
%     for j=2:length(Hr)
%         if Hr(j)<50 % on a choisi le seuil 20 pour determiner les petites region
%             Lr(Lr==j-1)=0; %suppression des petite regions -- L==j-1 récuprère l etiquete j-1
%         else
%             [yr,xr]=find(Lr==j-1);
%             vectr=[vectr; [mean(xr),mean(yr)]]; %calcul du centre de gravité des région
%                 % c est une concaténation vect=[vect; [mean(x),mean(y)]];
%         end 
%     end
%     
% %     figure(9),hold on, imagesc(Lr),colormap,title('image finale avec le centre de gravité du rouge');
% %     plot(vectr(:,1),vectr(:,2), '+r'); % trajectoire centre de gravité avec des carré rouge
%   
% % calcul du centre de gravité des carré vert
% Dg = I(:,:,1)< 15 &  I(:,:,2)>140 & I(:,:,3)>120 & I(:,:,3)<220;
% Dg=double(Dg);
% [Lg,numg]=bwlabel(Dg,4); %etiquetage région 4 connexité , L contient les regions  
% %numgegiong=numg; % affichage nombre de région
% Hg=hist(Lg(:),[0:numg]); 
%     for j=2:length(Hg)
%         if Hg(j)<51 % on a choisi le seuil 20 pour determiner les petites region
%             Lg(Lg==j-1)=0; %suppression des petite regions -- L==j-1 récuprère l etiquete j-1
%         else
%             [yg,xg]=find(Lg==j-1);
%             vectg=[vectg; [mean(xg),mean(yg)]]; %calcul du centre de gravité des région
%                 % c est une concaténation vect=[vect; [mean(x),mean(y)]];
%         end 
%     end
%     
% %     figure(10),hold on, imagesc(Lg),colormap,title('image finale avec le centre de gravité');
% %     plot(vectg(:,1),vectg(:,2), '+g'); % trajectoire centre de gravité avec des carré rouge
% 
% % caclul du centre de gravité du carré jaune    
%     
% Dy = I(:,:,1)> 190 &  I(:,:,2)> 190  & I(:,:,3)>100 & I(:,:,3)<210;
% Dy=double(Dy);
% [Ly,numy]=bwlabel(Dy,4); %etiquetage région 4 connexité , L contient les regions  
% %numgegiong=numy; % affichage nombre de région
% Hy=hist(Ly(:),[0:numy]); 
%     for j=2:length(Hy)
%         if Hy(j)<32 % on a choisi le seuil 20 pour determiner les petites region
%             Ly(Ly==j-1)=0; %suppression des petite regions -- L==j-1 récuprère l etiquete j-1
%         else
%             [yy,xy]=find(Ly==j-1);
%             vecty=[vecty; [mean(xy),mean(yy)]]; %calcul du centre de gravité des région
%                 % c est une concaténation vect=[vect; [mean(x),mean(y)]];
%         end 
%     end
%     
% %     figure(11),hold on, imagesc(Ly),colormap,title('image finale avec le centre de gravité');
% %     plot(vecty(:,1),vecty(:,2), '+y'); % trajectoire centre de gravité avec des carré rouge    
%     
% end

vecty = vecty(1:435,1:2);
vectb = vectb(1:435,1:2);

 figure(11),hold on,colormap,title('image finale avec le centre de gravité jaune');
 plot(vecty(:,1),vecty(:,2), '+r');
 
 figure(12),hold on,colormap,title('image finale avec le centre de gravité rouge');
 plot(vectr(:,1),vectr(:,2), '*r');
 
  figure(13),hold on,colormap,title('image finale avec le centre de gravité blue');
 plot(vectb(:,1),vectb(:,2), '+b');
 
  figure(14),hold on,colormap,title('image finale avec le centre de gravité vert');
 plot(vectg(:,1),vectg(:,2), '+g');
 
K = cameraParams.Intrinsics.IntrinsicMatrix'


% point World
pointr = [0 0 1]
pointrn = [0 0  0 1]
pointy = [187 0 1]
pointb = [0 109 1]
pointg = [186 108 1]
M0_homo = [0 0 0 1;187 0 0 1;0 109 0 1;186 108 0 1]

% % point images
% pro_point2 = cameraParams.ReprojectedPoints(2,:,1)
% pro_point24 = cameraParams.ReprojectedPoints(24,:,1)
% pro_point38 = cameraParams.ReprojectedPoints(38,:,1)
% pro_point10 = cameraParams.ReprojectedPoints(10,:,1)

HH = zeros(3, 4, 10);

for i= 1:10
% ensemble de point world
M0 = [pointr;pointy;pointb;pointg]
% ensemble de point image
M1 = [[vectr(i,1) vectr(i,2)] ;[vecty(i,1) vecty(i,2)];[vectb(i,1) vectb(i,2)];[vectg(i,1) vectg(i,2)]]

% calcul de l'homographie 9*1

% faut mettre M1 les point images et M0 les point world
H1= myDLT(M1,M0);

% matrice d'homographie 3*3
h1(1,1:3) = H1(1:3);
h1(2,1:3) = H1(4:6);
h1(3,1:3) = H1(7:9);




k_1H = inv(K)*h1;
alpha_r1 = k_1H(:,1);
alpha_r2 = k_1H(:,2);
alpha_t  = k_1H(:,3);

alpha2_r3 = cross(alpha_r1,alpha_r2);

alpha4 = det([alpha_r1 alpha_r2 alpha2_r3]);

alpha_ = -(alpha4)^(1/4);
alpha = (alpha4)^(1/4);
% calcul des vecteurs de la matrice de rotation

r1_ = alpha_r1 /alpha_;
r2_ = alpha_r2/alpha_;
r3_ = alpha2_r3 / (alpha_)^2;


r1 = alpha_r1 /alpha;
r2 = alpha_r2/alpha;
r3 = alpha2_r3 / (alpha)^2;

%format long


% matrice de rotation calculée
% R orienté suivant -z
R_ = [r1_ r2_ r3_]
% R orienté suivant z
R = [r1 r2 r3]

% matrice de rotation donnée par matlab
Rm =  cameraParams.RotationMatrices(:,:,1)'
% vecteur de translation calculé

% vecteur de translation
% suivant -z
t_ = alpha_t /alpha_
% suivant z
t = alpha_t /alpha

% vecteur de translation donné par matlab
tm = cameraParams.TranslationVectors(1,:)'


P = K*[R t]
PP(1:3,1:4,i) = P

Pm = K*[Rm tm]


end

pointrn = [0 0  50 1]
pointyn = [187 0 50 1]
pointbn = [0 109 50 1]
pointgn = [186 108 50 1]


pointr_cal=PP(:,:,1)*pointrn'

pointr_cal_n(1) = pointr_cal(1)/pointr_cal(3)
pointr_cal_n(2) = pointr_cal(2)/pointr_cal(3)

pointy_cal=PP(:,:,1)*pointyn'

pointy_cal_n(1) = pointy_cal(1)/pointy_cal(3)
pointy_cal_n(2) = pointy_cal(2)/pointy_cal(3)

pointb_cal=PP(:,:,1)*pointbn'

pointb_cal_n(1) = pointb_cal(1)/pointb_cal(3)
pointb_cal_n(2) = pointb_cal(2)/pointb_cal(3)


pointg_cal=PP(:,:,1)*pointgn'

pointg_cal_n(1) = pointg_cal(1)/pointg_cal(3)
pointg_cal_n(2) = pointg_cal(2)/pointg_cal(3)

I11 = imread('myVideoFrames2/000001.jpg');
figure
imagesc(I11)

hold on 
plot(pointr_cal_n(1),pointr_cal_n(2),'*r')
plot(pointy_cal_n(1),pointy_cal_n(2),'*r')
plot(pointb_cal_n(1),pointb_cal_n(2),'*r')
plot(pointg_cal_n(1),pointg_cal_n(2),'*r')

%point_img_calcul= zeros(
% point_img_calcul = PP(:,:,1)*M0_homo
% for i = 1:3
%     point_img_calcul2(i,1)= point_img_calcul(i,1)/