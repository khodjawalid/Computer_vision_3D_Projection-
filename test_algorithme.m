% Lire l'image
img = imread('AAA.jpg'); % Remplacez par le nom de votre image

% Convertir l'image en espace de couleur HSV
hsvImg = rgb2hsv(img);

% Définir les intervalles de la couleur rouge en HSV
hueMin = 0;    % Teinte minimale pour le rouge
hueMax = 0.05; % Teinte maximale pour le rouge
saturationMin = 0.5; % Saturation minimale
saturationMax = 1;   % Saturation maximale
valueMin = 0.2;  % Valeur minimale
valueMax = 1;    % Valeur maximale

% Créer un masque pour extraire la couleur rouge
redMask = (hsvImg(:,:,1) >= hueMin & hsvImg(:,:,1) <= hueMax) & ...
          (hsvImg(:,:,2) >= saturationMin & hsvImg(:,:,2) <= saturationMax) & ...
          (hsvImg(:,:,3) >= valueMin & hsvImg(:,:,3) <= valueMax);

% Appliquer le masque pour isoler les zones rouges
redRegion = redMask; % C'est déjà une image binaire

% Trouver les contours dans la région rouge
[contours, ~] = bwboundaries(redRegion, 'noholes');

% Afficher l'image originale
imshow(img); hold on;

% Rechercher les contours des carrés
for k = 1:length(contours)
    % Extraire le contour
    boundary = contours{k};

    % Calculer les longueurs des côtés du contour
    sideLengths = zeros(4, 1);
    for i = 1:4
        p1 = boundary(i, :);
        p2 = boundary(mod(i, 4) + 1, :);
        sideLengths(i) = sqrt((p2(1) - p1(1))^2 + (p2(2) - p1(2))^2);
    end
    
    % Vérifier si les longueurs des côtés sont égales (indiquant un carré)
    if max(sideLengths) - min(sideLengths) < 5 % Tolérance pour les longueurs similaires
        % Afficher le contour du carré rouge
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2); % Carré en vert
        % Afficher les coins du carré
        plot(boundary(1:4,2), boundary(1:4,1), 'ro', 'MarkerSize', 10); % Coins en rouge
    end
end

hold off;