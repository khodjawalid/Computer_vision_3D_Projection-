clear all
close all
load('cameraParams.mat')
%% Section une : Test du bon fonctionnement de la DLT et la matrice intrinsèque :

% extraction de la matrice intrinsèque à l'aide de caméra calibrator
K = cameraParams.Intrinsics.IntrinsicMatrix';

% Calcul des matrices extrinsèques à l'aide des homographies calculé par la fonction DLT

% point World
point2 = [cameraParams.WorldPoints(2,:) 1];
point24 = [cameraParams.WorldPoints(24,:) 1];
point38 = [cameraParams.WorldPoints(38,:) 1];
point10 = [cameraParams.WorldPoints(10,:) 1];

% point images correspondant
pro_point2 = cameraParams.ReprojectedPoints(2,:,1);
pro_point24 = cameraParams.ReprojectedPoints(24,:,1);
pro_point38 = cameraParams.ReprojectedPoints(38,:,1);
pro_point10 = cameraParams.ReprojectedPoints(10,:,1);

% ensemble de point world
M0 = [point2;point24;point38;point10];
% ensemble de point image
M1 = [pro_point2;pro_point24;pro_point38;pro_point10];

% calcul de l'homographie 9*1

% il faut mettre M1 les point images et M0 les point world
H1= myDLT(M1,M0);

% matrice d'homographie 3*3
h1(1,1:3) = H1(1:3);
h1(2,1:3) = H1(4:6);
h1(3,1:3) = H1(7:9);
h1;

% calcul des 2 matrices de projection possible
k_1H = inv(K)*h1;
alpha_r1 = k_1H(:,1);
alpha_r2 = k_1H(:,2);
alpha_t  = k_1H(:,3);


alpha2_r3 = cross(alpha_r1,alpha_r2);

% valeur du coefficient Alpha^4
alpha4 = det([alpha_r1 alpha_r2 alpha2_r3]);

% les deux valeur possible de alpha
alpha_ = -(alpha4)^(1/4);
alpha = (alpha4)^(1/4);

% calcul des vecteurs de la matrice de rotation

% avec le alpha_ négatif
r1_ = alpha_r1 /alpha_;
r2_ = alpha_r2/alpha_;
r3_ = alpha2_r3 / (alpha_)^2;

% avec le alpha_ positif
r1 = alpha_r1 /alpha;
r2 = alpha_r2/alpha;
r3 = alpha2_r3 / (alpha)^2;

%format long

% matrice de rotation calculée
% R orienté suivant -z
R_ = [r1_ r2_ r3_];
% R orienté suivant z
R = [r1 r2 r3];

% matrice de rotation donnée par matlab
Rm =  cameraParams.RotationMatrices(:,:,1)';


% vecteur de translation calculé
% suivant -z
t_ = alpha_t /alpha_;
% suivant z
t = alpha_t /alpha;

% vecteur de translation donné par matlab
tm = cameraParams.TranslationVectors(1,:)';

% les deux matrices de projections 
P = K*[R t];
P = K*[R_ t_];

% la matrice de projection donnée par matlab
Pm = K*[Rm tm];


% test et comparaison des résultats :
format short
N = P*[cameraParams.WorldPoints(2,:) 0 1; cameraParams.WorldPoints(24,:) 0 1;cameraParams.WorldPoints(38,:) 0 1;cameraParams.WorldPoints(10,:) 0 1]';

p1 = N(:,1);
disp("resultat calculée point 1")
p1n = p1(1:2)/p1(3);
disp("point image 1 donné par matlab")
pro_point2'; %%MATLAB


p2 = N(:,2);
disp("resultat calculée point 2")
p2n = p2(1:2)/p2(3);
disp("point image 2 donné par matlab")
pro_point24';

p3 = N(:,3);
disp("resultat calculée point 3")
p3n = p3(1:2)/p3(3);
disp("point image 3 donné par matlab")
pro_point38';

p4 = N(:,4);
disp("resultat calculée point 4")
p4n = p4(1:2)/p4(3);
disp("point image 4 donné par matlab")
pro_point10';

disp("les points images calculé par notre matrice de projection de  et ceux données par matlab sont les mêmes") 
%}

%% Section deux : génération d'un objet 3D dans la vidéo ( on génère une petite maison )



% récuperation des images de la video avec le fond noir
 %F = extractFrames("projet_fond_noir.mp4")
 % récupérées



% déclaration des vectuers qui vont sauvegarder la suite des centres de gravités
% de chacun des carrés

vectb = [];  % blue
vectr = [];  % red
vectg = [];  % green
vecty = [];  % yellow

%% calcul des centre de gravités des carrés dans toutes les images 
for n=001:435
I = imread(sprintf('myVideoFrames2/000%03d.jpg',n));

%
%----------------------------------------------------------------------------------------------------%  
% caclul du centre de gravité du carré bleu
Db = I(:,:,1)<40 &  I(:,:,2)>0 & I(:,:,2) < 130 & I(:,:,3) > 120; % application du filtre pour ne garder que le carré bleu
Db=double(Db);
[Lb,numb]=bwlabel(Db,4); % etiquetage région 4 connexité , Lb contient les regions bleu 

Hb=hist(Lb(:),[0:numb]);   % Histogramme de la photo résultante
    for j=2:length(Hb)
        if Hb(j)<110 % on a choisi le seuil 20 pour determiner les petites regions
            Lb(Lb==j-1)=0; % suppression des petites regions, Lb==j-1 récuprère l'etiquete j-1
        else    % si la région est asssez grande on en prend compte dans le calcul du centre d e gravité
            [yb,xb]=find(Lb==j-1);
            vectb=[vectb; [mean(xb),mean(yb)]]; %calcul du centre de gravité des région et ajout au vecteur 
               
        end 
    end
    
%     figure(8),hold on, imagesc(Lb),colormap,title('image finale avec le centre de gravité');
%     plot(vectb(:,1),vectb(:,2), '+b'); % trajectoire centre de gravité avec des carré bleu


%----------------------------------------------------------------------------------------------------%
%caclul du centre de gravité du carré rouge
    
Dr = I(:,:,1)>180 & I(:,:,2)>25 & I(:,:,2) < 150 &  I(:,:,3) > 55  & I(:,:,3) < 170;
Dr=double(Dr);
[Lr,numr]=bwlabel(Dr,4); % etiquetage région 4 connexité , Lr contient les regions  

Hr=hist(Lr(:),[0:numr]); 
    for j=2:length(Hr)
        if Hr(j)<50 % on a choisi le seuil 50 pour determiner les petites regions
            Lr(Lr==j-1)=0; % suppression des petite regions, Lr==j-1 récuprère l'etiquete j-1
        else
            [yr,xr]=find(Lr==j-1);
            vectr=[vectr; [mean(xr),mean(yr)]]; % calcul du centre de gravité des région et ajout au vecteur 
        end 
    end
    
%     figure(9),hold on, imagesc(Lr),colormap,title('image finale avec le centre de gravité du rouge');
%     plot(vectr(:,1),vectr(:,2), '+r'); % trajectoire centre de gravité avec des carré rouge

 %
%----------------------------------------------------------------------------------------------------%
% calcul du centre de gravité des carré vert

Dg = I(:,:,1)<15 &  I(:,:,2)>140 & I(:,:,3)>120 & I(:,:,3)<220;
Dg=double(Dg);
[Lg,numg]=bwlabel(Dg,4); % etiquetage région 4 connexité , L contient les regions  

Hg=hist(Lg(:),[0:numg]);
    for j=2:length(Hg)
        if Hg(j)<40 % on a choisi le seuil 50 pour determiner les petites regions
            Lg(Lg==j-1)=0; % suppression des petite regions; Lg==j-1 récuprère l'etiquete j-1
        else
            [yg,xg]=find(Lg==j-1);
            vectg=[vectg; [mean(xg),mean(yg)]]; % calcul du centre de gravité des région et ajout au vecteur 
        end
    end
    
%     figure(10),hold on, imagesc(Lg),colormap,title('image finale avec le centre de gravité');
%     plot(vectg(:,1),vectg(:,2), '+g'); % trajectoire centre de gravité avec des carré rouge

%----------------------------------------------------------------------------------------------------%
% caclul du centre de gravité du carré jaune    
  
Dy = I(:,:,1)> 190 &  I(:,:,2)> 190  & I(:,:,3)>100 & I(:,:,3)<210;
Dy=double(Dy);
[Ly,numy]=bwlabel(Dy,4); % etiquetage région 4 connexité , Ly contient les regions  
Hy=hist(Ly(:),[0:numy]); 
    for j=2:length(Hy)
        if Hy(j)<32 % on a choisi le seuil 32 pour determiner les petites regions
            Ly(Ly==j-1)=0; %suppression des petite regions, Ly==j-1 récuprère l'etiquete j-1
        else
            [yy,xy]=find(Ly==j-1);
            vecty=[vecty; [mean(xy),mean(yy)]]; %calcul du centre de gravité des région et ajout au vecteur   
        end 
    end
    
%     figure(11),hold on, imagesc(Ly),colormap,title('image finale avec le centre de gravité');
%     plot(vecty(:,1),vecty(:,2), '+y'); % trajectoire centre de gravité avec des carré rouge    
    
 end

% je ne  garde que les centres de gravité dans 435 images car ici on a 2 ou 3
% Centre de gravité en plus du à l'erreur de detection des carré
vecty = vecty(1:435,1:2);
vectb = vectb(1:435,1:2);
vectg = vectg(1:435,1:2);
vectr = vectr(1:435,1:2);


% tracé des trajectoires des centres de gravités des 4 carrés sur 4 figures
% différentes

 figure(11),hold on,colormap,title('image finale avec le centre de gravité jaune');
 plot(vecty(:,1),vecty(:,2), '+r');
 
 figure(12),hold on,colormap,title('image finale avec le centre de gravité rouge');
 plot(vectr(:,1),vectr(:,2), '*r');
 
  figure(13),hold on,colormap,title('image finale avec le centre de gravité blue');
 plot(vectb(:,1),vectb(:,2), '+b');
  
  figure(14),hold on,colormap,title('image finale avec le centre de gravité vert');
 plot(vectg(:,1),vectg(:,2), '+g');
 
 
 
 
 
%}

 

%% calcul des matrices de projections pour mes différentes images

%
%point World
pointr = [0 0 1]
pointy = [187 0 1]
pointb = [0 109 1]
pointg = [186 108 1]
M0_homo = [0 0 0 1;187 0 0 1;0 109 0 1;186 108 0 1]


% Matrice à 3 dimension pour stocker toutes les matrices de projection
PP = zeros(3, 4, 435);
K = cameraParams.Intrinsics.IntrinsicMatrix'
for i= 1:435
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

% vecteur de translation calculé

% vecteur de translation
% suivant -z
t_ = alpha_t /alpha_
% suivant z
t = alpha_t /alpha

% les matrices de projection
P = K*[R t]
P_ = K*[R_ t_]
% choix de la matrice de projection à utiliser pour la génération de la
% maison
if t(3) < 0
    PP(1:3,1:4,i) = P
else 
    PP(1:3,1:4,i) = P_
end

end

%
%% génération de la maison
%

% les vecteur qui vont stocker les coordonnées des point de la maison
pointrn = zeros(10,4)
pointyn = zeros(10,4)
pointbn = zeros(10,4)
pointgn = zeros(10,4)

pointry = zeros(37,4)
pointrb = zeros(21,4)
pointgy = zeros(21,4)
pointbg = zeros(37,4)

pointrg = zeros(21,4)

pointby = zeros(21,4)

% les lignes verticales de la maison dans le repere monde
for i = 1:11
    
pointrn(i,:) = [0 0 (i-1)*5  1]
pointyn(i,:) = [187 0 (i-1)*5 1]
pointbn(i,:) = [0 109 (i-1)*5 1]
pointgn(i,:) = [186 108 (i-1)*5 1]
end

%les lignes horizontale dans le repere monde
% entre le rouge et le jaune
for i= 1:38   
pointry(i,:) = [2+(i-1)*5 0 50  1]
end 

% entre le rouge et le bleu
for i= 1:22  
pointrb(i,:) = [0 4+(i-1)*5 50  1]

end 
% entre le vert et le jaune
for i= 1:22  
pointgy(i,:) = [187  3+(i-1)*5 50 1]
end 
% entre le bleu et le vert
for i= 1:38   
pointbg(i,:) = [1+(i-1)*5 109 50 1]
end
%entre le rouge et le vert pour le toit de la maison
for i= 1:22
    if i <12
        pointrg(i,:) = [round((186*5*(i-1))/108) 3+(i-1)*5 50+2*i  1]
    else
        pointrg(i,:) = [round((186*5*(i-1))/108) 3+(i-1)*5 72-2*(i-11)  1]
    end
end
for i= 1:22
    if i <12
        pointby(i,:) = [round((-186*5*(i-1))/108)+186 3+(i-1)*5 50+2*i  1]
    else
        pointby(i,:) = [round((-186*5*(i-1))/108)+186 3+(i-1)*5 72-2*(i-11)  1]
    end
end
 
for i = 1:150
I11 = imread(sprintf('myVideoFrames2/000%03d.jpg',i));
figure(1)
imagesc(I11)
for j = 1:10
%calcul des projections ( les points images)
pointr_cal=PP(:,:,i)*pointrn(j,:)';

pointr_cal_n(1) = pointr_cal(1)/pointr_cal(3);
pointr_cal_n(2) = pointr_cal(2)/pointr_cal(3);

pointy_cal=PP(:,:,i)*pointyn(j,:)';

pointy_cal_n(1) = pointy_cal(1)/pointy_cal(3);
pointy_cal_n(2) = pointy_cal(2)/pointy_cal(3);

pointb_cal=PP(:,:,i)*pointbn(j,:)'

pointb_cal_n(1) = pointb_cal(1)/pointb_cal(3);
pointb_cal_n(2) = pointb_cal(2)/pointb_cal(3);


pointg_cal=PP(:,:,i)*pointgn(j,:)';

pointg_cal_n(1) = pointg_cal(1)/pointg_cal(3);
pointg_cal_n(2) = pointg_cal(2)/pointg_cal(3);



hold on 
plot(pointr_cal_n(1),pointr_cal_n(2),'*r')
plot(pointy_cal_n(1),pointy_cal_n(2),'*r')
plot(pointb_cal_n(1),pointb_cal_n(2),'*r')
plot(pointg_cal_n(1),pointg_cal_n(2),'*r')
end

for j = 1:38
pointr_cal=PP(:,:,i)*pointry(j,:)';

pointr_cal_n(1) = pointr_cal(1)/pointr_cal(3);
pointr_cal_n(2) = pointr_cal(2)/pointr_cal(3);

pointy_cal=PP(:,:,i)*pointbg(j,:)';

pointy_cal_n(1) = pointy_cal(1)/pointy_cal(3);
pointy_cal_n(2) = pointy_cal(2)/pointy_cal(3);
if j<22
    pointb_cal=PP(:,:,i)*pointrb(j,:)';

    pointb_cal_n(1) = pointb_cal(1)/pointb_cal(3);
    pointb_cal_n(2) = pointb_cal(2)/pointb_cal(3);

    pointg_cal=PP(:,:,i)*pointgy(j,:)';

    pointg_cal_n(1) = pointg_cal(1)/pointg_cal(3);
    pointg_cal_n(2) = pointg_cal(2)/pointg_cal(3);
    %ligne toit
    pointrg_cal=PP(:,:,i)*pointrg(j,:)';

    pointrg_cal_n(1) = pointrg_cal(1)/pointrg_cal(3);
    pointrg_cal_n(2) = pointrg_cal(2)/pointrg_cal(3);
    %ligne toit 2 
    pointby_cal=PP(:,:,i)*pointby(j,:)';

    pointby_cal_n(1) = pointby_cal(1)/pointby_cal(3);
    pointby_cal_n(2) = pointby_cal(2)/pointby_cal(3);

end
hold on 
plot(pointr_cal_n(1),pointr_cal_n(2),'*r')
plot(pointy_cal_n(1),pointy_cal_n(2),'*r')
if j<22
    plot(pointb_cal_n(1),pointb_cal_n(2),'*r')
    plot(pointg_cal_n(1),pointg_cal_n(2),'*r')
    %plot du toit
    plot(pointrg_cal_n(1),pointrg_cal_n(2),'*r')
    plot(pointby_cal_n(1),pointby_cal_n(2),'*r')
    

end
end
end

%}

