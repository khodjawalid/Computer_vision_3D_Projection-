%% Calibration et calcul de la matrice H 

%Mesurre des coordonnées des points en pixels et en mm

p1 = [0 0;0 5;8 0;8 5;1 4;3 2;5 4;7 1;3 5;3 0]*21; % [mm]
p2 = [405 727;409 908;761 704;818 873;458 862;550 780;659 846;726 735;567 892;544 718]; % [px]

%Extractionn du nombre de points
n = size(p1,1);

%Initialisation de la matrice H 
A = zeros(2*n, 9);

%Remplissage de la matrice A 
for i = 1:n
    x = p1(i,1);
    y = p1(i,2);
    xp = p2(i,1);
    yp = p2(i,2);

    A(2*i-1,:) = [0 0 0 -x -y -1 x*yp y*yp yp];
    A(2*i,:) = [x y 1 0 0 0 -x*xp -y*xp -xp];
    
end

%Utilisation de la méthode DLT 
[U,S,V] = svd(A);
h = V(:,end);
H = reshape(h,[3 3])';
 
% Affichage de l'image
img = imread('Points2.jpeg'); 

figure(1); % Create a new figure window
imshow(img); % Display the image
hold on;

%Méthode de hadriel sans stockage de doonnées 
%Calcul des coordonnées des nouveaux points et les afficher 
% for i = 1:n
%     xr = p1(i,:);
%     xpt = H*[p1(i,:) 1]'*21;
%     xpt = xpt/xpt(3);
%     plot(xpt(1),xpt(2),'ro')
% end

%Ma méthode avec sauvegarde de données 
%Calcul des coordonnées des nouveaux points et les afficher 

xpt1=zeros(n,3);
for i = 1:n
    xpt1(i,:) = H*[p1(i,:) 1]'*21;
    xpt1(i,:) = xpt1(i,:)/xpt1(i,3);
end

x = xpt1(:, 1); % Extraire les coordonnées X
y = xpt1(:, 2); % Extraire les coordonnées Y

plot(x, y, 'bo'); % 'b' pour bleu et 'o' pour des marqueurs circulaires
xlabel('X');
ylabel('Y');
title("Affichage de l'image avec les points respectifs");


%% Affichage d'un carrée en partant des coordonnées réelles

% Coordonnées des extrémités du carrée 
x_min = 1;
y_min = 1;
x_max = 6;
y_max = 4;

% Définir les sommets du carré
xcar = [x_min, x_max, x_max, x_min, x_min]*21; % [mm]
ycar = [y_min, y_min, y_max, y_max, y_min]*21; %[mm]

n = length(xcar);
ptcar=zeros(n,3);

%Calculer les données des extrimitées en pixels 
for i = 1:n
    ptcar(i,:) = H*[xcar(i) ycar(i) 1]'*21;
    ptcar(i,:) = ptcar(i,:)/ptcar(i,3);
end

xcarpix = ptcar(:, 1); % Extraire les coordonnées X
ycarpix = ptcar(:, 2); % Extraire les coordonnées Y

%Afficher le carrée 
figure(2);
imshow(img); % Display the image
hold on;

% Tracer le carré
plot(xcarpix, ycarpix,'-o' ,'LineWidth', 2, 'MarkerSize', 8); % '-o' pour les lignes avec marqueurs
title('Carré à partir des coordonnées des extrémités');

%% Calcul de la mattrice H 

k = [1.1546e3, 0, 0.5945e3; 0, 1.1537e3, 0.8078; 0, 0, 0.0010e3];

kinv = inv(k);

RT = kinv * H;

r1 = RT(:,1);
r2 = RT(:,2);
t= RT(:,3);

r1 = r1/norm(r1);
r2 = r2/norm(r2);

r3 = cross(r1,r2);

R = [r1,r2,r3];

trans_h = [R, t];

P = k*trans_h



