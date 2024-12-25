%% Calibration et calcul de la matrice H 


%focale 0.6


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

    %A(2*i-1,:) = [0 0 0 -x -y -1 x*yp y*yp yp];
    %A(2*i,:) = [x y 1 0 0 0 -x*xp -y*xp -xp];

    A(2*i,:) = [0 0 0 x y 1 -x*yp -y*yp -yp];
    A(2*i-1,:) = [x y 1 0 0 0 -x*xp -y*xp -xp];
    
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

%% Calcul de la mattrice P 

k = [1.1546e3  , 0 ,     0.5945e3; 
      0    , 1.1537e3  , 0.8078e3;
      0    ,     0     , 0.0010e3];

kinv = inv(k);

Hn=H/H(end,end);

RT = kinv * Hn;

r1 = RT(:,1);
r2 = RT(:,2);
t= RT(:,3);
r3 = cross(r1,r2);

R = [r1,r2,r3];

alpha = nthroot(det(R),4);

r1 = r1/alpha;
r2 = r2/alpha;
t=t/alpha;
r3=r3/alpha^2;

R = [r1,r2,r3];

%On doit verifier la rotation et la translation positif ou négative 
Rn =[-R(:,1),-R(:,2), R(:,3)];
tn=-t
trans_h = [Rn, tn];

P = alpha*k*trans_h;

%% Test de la fonction P 

mreel1 = [1;1;2]*21; %[mm]
mreel2 = [1;1;0]*21;

mpxl1 = P * [mreel1;1];
mpxl2 = P * [mreel2;1];

mpxl1 = mpxl1/mpxl1(end);
mpxl2 = mpxl2/mpxl2(end);

%figure(3);

% imshow(img); % Display the image
%hold on;

% Tracer le carré
%plot(mpxl1(1), mpxl1(2),'bo'); % '-o' pour les lignes avec marqueurs
%plot(mpxl2(1), mpxl2(2),'bo'); % '-o' pour les lignes avec marqueurs


% Définir les sommets du cube dans le repère réel (en mm)
cube_real = [
    0, 0, 0;  % S1
    0, 0, 5; % S2
    5, 0, 0; % S3
    5, 0, 5; % S4
    0, 5, 0; % S5
    0, 5, 5; % S6
    5, 5, 0; % S7
    5, 5, 5; % S8
] * 21; % Échelle en mm (exemple)

% Ajouter la coordonnée homogène (1) à chaque sommet
cube_real_h = [cube_real, ones(size(cube_real, 1), 1)];

% Projeter les sommets dans le repère pixels
cube_pixel_h = (P * cube_real_h')'; % Projeter avec la matrice P
cube_pixel_h = cube_pixel_h ./ cube_pixel_h(:, 3); % Normaliser (diviser par la 3e coordonnée)

% Extraire les coordonnées X et Y en pixels
cube_pixel = cube_pixel_h(:, 1:2);

% Définir les arêtes du cube (indices des sommets à relier)
edges = [
    1, 2; 1, 3; 1, 5; % Arêtes depuis S1
    2, 4; 2, 6;       % Arêtes depuis S2
    3, 4; 3, 7;       % Arêtes depuis S3
    4, 8;             % Arête depuis S4
    5, 6; 5, 7;       % Arêtes depuis S5
    6, 8;             % Arête depuis S6
    7, 8;             % Arête depuis S7
];

% Affichage de l'image et du cube
figure;
imshow(img);
hold on;

% Tracer les arêtes du cube
for i = 1:size(edges, 1)
    % Récupérer les indices des sommets
    idx1 = edges(i, 1);
    idx2 = edges(i, 2);
    
    % Récupérer les coordonnées des deux sommets
    pt1 = cube_pixel(idx1, :);
    pt2 = cube_pixel(idx2, :);
    
    % Tracer une ligne entre les deux sommets
    plot([pt1(1), pt2(1)], [pt1(2), pt2(2)], 'r-', 'LineWidth', 2);
end

% Afficher les sommets du cube
plot(cube_pixel(:, 1), cube_pixel(:, 2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

hold off;




%Blinder communiity


