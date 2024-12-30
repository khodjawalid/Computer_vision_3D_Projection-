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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calcul de la mattrice P en 3d
    
    %Matrice intinsèque 
    k = [1.1546e3  , 0 ,     0.5945e3; 
          0    , 1.1537e3  , 0.8078e3;
          0    ,     0     , 0.0010e3];

    %Calcul de RT 
    kinv = inv(k);
    Hn=H/H(end,end);
    RT = kinv * Hn;
    
    % Extraire rotation et translation
    r1 = RT(:,1);
    r2 = RT(:,2);
    t= RT(:,3);
    r3 = cross(r1,r2);
    R = [r1,r2,r3];

    %Coefficient de normalisation 
    alpha = nthroot(det(R),4);
    
    %Normalisation
    r1 = r1/alpha;
    r2 = r2/alpha;
    t=t/alpha;
    r3=r3/alpha^2;    
    R = [r1,r2,r3];
    
    %On doit verifier la rotation et la translation positif ou négative 
    Rn =[-R(:,1),-R(:,2), -R(:,3)];
    tn=-t;
    
    %Nouvelle matrice de transformation homogène 
    trans_h = [Rn, tn];
    
    % Construire la matrice de projection
    P = alpha*k*trans_h; 
    
    eiffel = vertices*80 ; % Échelle en mm 
    eiffel(:, 1:2) =eiffel(:, 1:2) +50; %Centrer la tour eiffel au milieu 
    
    % Ajouter la coordonnée homogène (1) à chaque sommet
    eiffel_h = [eiffel, ones(size(eiffel, 1), 1)];
    
    % Calculer les couleurs basées sur la valeur de x en coordonnées réelles
    % Calculer les couleurs basées sur la valeur de y en coordonnées réelles
    colors = zeros(size(eiffel, 1), 3); % Initialiser une matrice pour les couleurs
    for i = 1:size(eiffel, 1)
        if eiffel(i, 2) < 40
            colors(i, :) = [0, 0, 1]; % Bleu
        elseif eiffel(i, 2) < 65
            colors(i, :) = [1, 1, 1]; % Blanc
        else
            colors(i, :) = [1, 0, 0]; % Rouge
        end
    end

    % Projeter les sommets dans le repère pixels
    eiffel_pixel_h = (P * eiffel_h')'; % Projeter avec la matrice P
    eiffel_pixel_h = eiffel_pixel_h ./ eiffel_pixel_h(:, 3); % Normaliser (diviser par la 3e coordonnée)
    
    % Extraire les coordonnées X et Y en pixels
    eiffel_pixel = eiffel_pixel_h(:, 1:2);
    
    % Affichage de l'image et du cube
    imshow(I);
    hold on;
    
    % Afficher les sommets du cube
    scatter(eiffel_pixel(:, 1), eiffel_pixel(:, 2), 5, colors, 'filled');
    
    hold off;
    
    pause(1/40)
end



