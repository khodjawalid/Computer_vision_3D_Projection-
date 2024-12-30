clear all;
clc;

%% Paramètres du filtre manuel pour la détection de couleur bleue
hueThresholdLow = 0.55; % Borne basse de la teinte
hueThresholdHigh = 0.75; % Borne haute de la teinte
saturationThreshold = 0.2; % Seuil minimum pour la saturation
valueThreshold = 0.2; % Seuil minimum pour la luminosité

% Initialisation des erreurs
errors = [];

% Boucle sur toutes les images (20 à 48)
for i = 20:48

    % Charger l'image
    I = imread(['Set3/set3_img (', num2str(i), ').jpg']);

    % Vérifier et corriger si nécessaire l'orientation de la photo
    info = imfinfo(['Set3/set3_img (', num2str(i), ').jpg']);
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

    %% Méthode 1 : Détection manuelle par segmentation HSV
    % Convertir en espace HSV
    hsvImage = rgb2hsv(I);
    hue = hsvImage(:, :, 1);        % Teinte
    saturation = hsvImage(:, :, 2); % Saturation
    value = hsvImage(:, :, 3);      % Valeur

    % Créer un masque binaire pour les pixels bleus
    blueMask = (hue >= hueThresholdLow) & (hue <= hueThresholdHigh) & ...
               (saturation >= saturationThreshold) & ...
               (value >= valueThreshold);

    % Nettoyer le masque avec des opérations morphologiques
    blueMask = imopen(blueMask, strel('square', 3)); % Suppression de bruit
    blueMask = imclose(blueMask, strel('square', 3)); % Remplir les trous

    % Calcul des régions et des centres de gravité
    stats = regionprops(blueMask, 'Centroid');
    manualCenters = cat(1, stats.Centroid);

    %% Méthode 2 : Détection automatisée avec des points d'intérêt
    grayImage = rgb2gray(I);
    points = detectSURFFeatures(grayImage);
    automatedCenters = points.Location;

    %% Comparaison des centres de gravité
    % Calcul des erreurs entre les centres détectés
    numPoints = min(size(manualCenters, 1), size(automatedCenters, 1));
    manualCenters = manualCenters(1:numPoints, :);
    automatedCenters = automatedCenters(1:numPoints, :);

    error = sqrt(sum((manualCenters - automatedCenters).^2, 2));
    errors = [errors; error];

    %% Affichage des résultats
    imshow(I);
    hold on;
    plot(manualCenters(:, 1), manualCenters(:, 2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(automatedCenters(:, 1), automatedCenters(:, 2), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    legend('Centres Manuels', 'Centres Automatisés');
    title(['Image ', num2str(i), ' - Comparaison des Centres']);
    hold off;
    pause(0.5);
end

%% Analyse finale des erreurs
meanError = mean(errors);
stdError = std(errors);

fprintf('Erreur moyenne : %.2f pixels\n', meanError);
fprintf('Écart-type des erreurs : %.2f pixels\n', stdError);