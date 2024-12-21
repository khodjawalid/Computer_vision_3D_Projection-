%% test
%Téléchargement de l'image

%Paramètres du filtre 
% Définir les seuils pour détecter le bleu (ajuster si nécessaire)
hueThresholdLow = 0.55; % Borne basse de la teinte
hueThresholdHigh = 0.75; % Borne haute de la teinte
saturationThreshold = 0.2; % Seuil minimum pour la saturation
valueThreshold = 0.2; % Seuil minimum pour la luminosité


for i=71:72 

    
    I = imread(['Frames2/IMG-20241221-WA00',num2str(i),'.jpg']);
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

    
    p_calcul = c
    p1 = [0 0; 0 5; 0 10; 5 0; 5 10; 10 0; 10 5; 10 10]*10; % [mm]10
    p_inter = zeros(8,2);


    if i == 71 
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
 
 
    close all;
    figure; % Create a new figure window
    imshow(I); % Display the image
    hold on;

    xpt1=zeros(n,3);
    for j = 1:n
        xpt1(j,:) = H*[p1(j,:) 1]'*21;
        xpt1(j,:) = xpt1(j,:)/xpt1(j,3);
    end

    x = xpt1(:, 1); % Extraire les coordonnées X
    y = xpt1(:, 2); % Extraire les coordonnées Y

    plot(x, y, 'ro'); % 'b' pour bleu et 'o' pour des marqueurs circulaires
    xlabel('X');
    ylabel('Y');
    title("Affichage de l'image avec les points respectifs");
    pause(0.5)

end
