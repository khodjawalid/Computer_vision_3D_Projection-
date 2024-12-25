%% test
%Téléchargement de l'image

%Paramètres du filtre 
% Définir les seuils pour détecter le bleu (ajuster si nécessaire)
hueThresholdLow = 0.55; % Borne basse de la teinte
hueThresholdHigh = 0.75; % Borne haute de la teinte
saturationThreshold = 0.2; % Seuil minimum pour la saturation
valueThreshold = 0.2; % Seuil minimum pour la luminosité


for i=1:48

    
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

    
    p_calcul = c;
    p1 = [0 0; 5 0; 0 5; 10 0; 0 10; 5 10; 10 5; 10 10]*10; % [mm]10
    p_inter = zeros(8,2);


    if i == 1 
        p_inter = p_calcul;  %pour la premières itération seulement 
    else 
        for j = 1:8 
            distances = sqrt(sum((p_sorted(j,:) - p_calcul).^2, 2));
            % Trouver l'indice du point le plus proche
            [~, idx] = min(distances);
            % Extraire le point le plus proche
            p_inter(j,:) = p_calcul(idx,:);
        end
    end 

    p_sorted = p_inter;


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

    %A(2*i-1,:) = [0 0 0 -x -y -1 x*yp y*yp yp];
    %A(2*i,:) = [x y 1 0 0 0 -x*xp -y*xp -xp];

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
    
    P = alpha*k*trans_h;

    % Test de la fonction P 

    mreel1 = [1;1;2]*21; %[mm]
    mreel2 = [1;1;0]*21;
    
    mpxl1 = P * [mreel1;1];
    mpxl2 = P * [mreel2;1];
    
    mpxl1 = mpxl1/mpxl1(end);
    mpxl2 = mpxl2/mpxl2(end);
    
    
    % Définir les sommets du cube dans le repère réel (en mm)
    cube_real = [
        0, 0, 0;  % S1
        0, 0, 5; % S2
        10, 0, 0; % S3
        10, 0, 5; % S4
        0, 10, 0; % S5
        0, 10, 5; % S6
        10, 10, 0; % S7
        10, 10, 5; % S8
    ] * 10; % Échelle en mm (exemple)
    
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
    
    imshow(I);
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
    pause(1/10)
    

end
