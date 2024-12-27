clear all
clc
%% Script de traitement d'image et projection 3D
% Ce script effectue les étapes suivantes :
% 1. Télécharge et corrige l'orientation des images.
% 2. Détecte les régions bleues dans chaque image.
% 3. Calcule une transformation homographique pour aligner les points.
% 4. Projette un cube 3D sur l'image à l'aide d'une matrice de projection.


% Paramètres du filtre pour la détection de couleur bleue
hueThresholdLow = 0.55; % Borne basse de la teinte
hueThresholdHigh = 0.75; % Borne haute de la teinte
saturationThreshold = 0.2; % Seuil minimum pour la saturation
valueThreshold = 0.2; % Seuil minimum pour la luminosité
figure()

% Charger le modèle 3D de la Tour Eiffel (format STL)
model = stlread('Eiffel.stl');
vertices = model.Points;
faces = model.ConnectivityList;


% Créer un objet VideoWriter
video = VideoWriter('animation', 'MPEG-4'); % 'animation' est le nom de la vidéo
open(video); % Ouvre le fichier vidéo pour écrire

% Boucle sur toutes les images (48 au total)
for i=20:48

    i

    % Charger l'image
    I = imread(['Set3/set3_img (',num2str(i),').jpg']);
    
    %Vérifier et corriger si nécessaire l'orientation de la photo 
    info = imfinfo(['Set3/set3_img (',num2str(i),').jpg']);
    if isfield(info, 'Orientation')
        switch info.Orientation
            case 3
                I = imrotate(I, 180); % Tourner de 180°
            case 6
                I = imrotate(I, -90); % Tourner de 90° antihoraire
            case 8
                I = imrotate(I, 90); % Tourner de 90° horaire
        end
    end

    % Convertir en espace HSV
    hsvImage = rgb2hsv(I);

    % Séparer les canaux HSV
    hue = hsvImage(:,:,1);        % Teinte
    saturation = hsvImage(:,:,2); % Saturation
    value = hsvImage(:,:,3);      % Valeur

    % Créer un masque binaire pour les pixels bleus
    blueMask = (hue >= hueThresholdLow) & (hue <= hueThresholdHigh) & ...
           (saturation >= saturationThreshold) & ...
           (value >= valueThreshold);

    % Optionnel : Nettoyer le masque avec des opérations morphologiques
    blueMask = imopen(blueMask, strel('square', 3)); % Suppression de bruit
    blueMask = imclose(blueMask, strel('square', 3)); % Remplir les trous
    
    % Créer l'image en noir et blanc 
    outputImage = double(blueMask) * 255; % Convertir le masque en image 8 bits

    % etiquetage région 4 connexité , Lb contient les regions bleu
    [Ir,numb]=bwlabel(outputImage); 
    s=zeros(numb,1); %stocker la taille 
    c=zeros(numb,2); %pour stocker les coordonnées x,y

    for j=1:numb 
        [y,x ]=find(Ir==j);
        c(j,:)=[mean(x),mean(y)];
        s(j)=length(x);
    end

    % Points calculés à partir des régions détectées
    p_calcul = c;
    p1 = [0 0; 5 0; 10 0; 0 5; 10 5; 0 10; 5 10; 10 10]*10; % [mm]10
    p_inter = zeros(8,2);

    % Initialisation des correspondances
    if i == 20 
        p_inter = p_calcul;  % Première image
    else 
        for j = 1:8 
            % Associer les points détectés aux points réels par distance minimale
            distances = sqrt(sum((p_sorted(j,:) - p_calcul).^2, 2));
            % Trouver l'indice du point le plus proche
            [~, idx] = min(distances);
            % Extraire le point le plus proche
            p_inter(j,:) = p_calcul(idx,:);
        end
    end 

    p_sorted = p_inter;

    % Calcul de la transformation homographique (DLT)
    p2 = p_sorted;

    %Extractionn du nombre de points
    n = size(p1,1);

    %Initialisation de la matrice H 
    A = zeros(2*n, 9);

    %Remplissage de la matrice A 
    for j = 1:n
        x = p1(j,1);
        y = p1(j,2);
        xp = p2(j,1);
        yp = p2(j,2);
        A(2*j,:) = [0 0 0 x y 1 -x*yp -y*yp -yp];
        A(2*j-1,:) = [x y 1 0 0 0 -x*xp -y*xp -xp];
    
    end

    %Utilisation de la méthode DLT 
    [U,S,V] = svd(A);
    h = V(:,end);
    H = reshape(h,[3 3])';
 
    %close all;
    % figure; % Create a new figure window
    % imshow(I); % Display the image
    % hold on;

    xpt1=zeros(n,3);
    for j = 1:n
        xpt1(j,:) = H*[p1(j,:) 1]'*21;
        xpt1(j,:) = xpt1(j,:)/xpt1(j,3);
    end

    x = xpt1(:, 1); % Extraire les coordonnées X
    y = xpt1(:, 2); % Extraire les coordonnées Y

    % plot(x, y, 'ro'); % 'b' pour bleu et 'o' pour des marqueurs circulaires
    % xlabel('X');
    % ylabel('Y');
    % title("Affichage de l'image avec les points respectifs");
    % pause(0.5)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calcul de la mattrice P 

    k = [1.1546e3  , 0 ,     0.5945e3; 
          0    , 1.1537e3  , 0.8078e3;
          0    ,     0     , 0.0010e3];
    kinv = inv(k);
    Hn=H/H(end,end);
    RT = kinv * Hn;
    
    % Extraire rotation et translation
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
    Rn =[-R(:,1),-R(:,2), -R(:,3)];
    tn=-t;
    
    trans_h = [Rn, tn];
    
    % Construire la matrice de projection
    P = alpha*k*trans_h;

    % Test de la fonction P 
    
   
    % Définir les sommets du cube dans le repère réel (en mm)

    
  
    cube_real = vertices*80 ; % Échelle en mm (exemple)
    cube_real(:, 1:2) =cube_real(:, 1:2) +50;
    
    % Ajouter la coordonnée homogène (1) à chaque sommet
    cube_real_h = [cube_real, ones(size(cube_real, 1), 1)];
    
    % Projeter les sommets dans le repère pixels
    cube_pixel_h = (P * cube_real_h')'; % Projeter avec la matrice P
    cube_pixel_h = cube_pixel_h ./ cube_pixel_h(:, 3); % Normaliser (diviser par la 3e coordonnée)
    
    % Extraire les coordonnées X et Y en pixels
    cube_pixel = cube_pixel_h(:, 1:2);
    
    % Affichage de l'image et du cube
    imshow(I);
    hold on;
    
    % Afficher les sommets du cube
    plot(cube_pixel(:, 1), cube_pixel(:, 2), 'o', 'MarkerSize', 2.5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    
    hold off;
    
    pause(1/40)
    
    

end



