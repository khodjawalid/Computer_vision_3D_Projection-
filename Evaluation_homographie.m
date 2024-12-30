clear all;
clc;

% Paramètres pour le filtre de détection de bleu
hueThresholdLow = 0.55; % Borne basse de la teinte
hueThresholdHigh = 0.75; % Borne haute de la teinte
saturationThreshold = 0.2; % Seuil minimum pour la saturation
valueThreshold = 0.2; % Seuil minimum pour la luminosité

% Points réels de référence pour la calibration
p1 = [0 0; 5 0; 10 0; 0 5; 10 5; 0 10; 5 10; 10 10] * 10; % [mm]

% Initialisation des erreurs
errors_manual = [];
errors_matlab = [];

for i = 20:48
    fprintf("Traitement de l'image",i,'\n');

    % Charger l'image
    I = imread(['Set3/set3_img (', num2str(i), ').jpg']);
    
    % Vérifier et corriger l'orientation si nécessaire
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

    % Convertir en espace HSV
    hsvImage = rgb2hsv(I);
    hue = hsvImage(:,:,1);        % Teinte
    saturation = hsvImage(:,:,2); % Saturation
    value = hsvImage(:,:,3);      % Valeur

    % Détecter les pixels bleus
    blueMask = (hue >= hueThresholdLow) & (hue <= hueThresholdHigh) & ...
               (saturation >= saturationThreshold) & ...
               (value >= valueThreshold);

    % Nettoyer le masque
    blueMask = imopen(blueMask, strel('square', 3));
    blueMask = imclose(blueMask, strel('square', 3));

    % Etiquetage des régions détectées
    [Ir, numb] = bwlabel(blueMask);
    c = zeros(numb, 2); % Stocker les coordonnées moyennes
    for j = 1:numb
        [y, x] = find(Ir == j);
        c(j, :) = [mean(x), mean(y)];
    end

    % Correspondances des points détectés
    if i == 20
        p_inter = c; % Initialisation
    else
        p_inter = zeros(size(p1));
        for j = 1:size(p1, 1)
            distances = sqrt(sum((c - p_sorted(j, :)).^2, 2));
            [~, idx] = min(distances);
            p_inter(j, :) = c(idx, :);
        end
    end
    p_sorted = p_inter;

    % Homographie manuelle (DLT)
    n = size(p1, 1);
    A = zeros(2 * n, 9);
    for j = 1:n
        x = p1(j, 1);
        y = p1(j, 2);
        xp = p_sorted(j, 1);
        yp = p_sorted(j, 2);
        A(2*j-1, :) = [x y 1 0 0 0 -x*xp -y*xp -xp];
        A(2*j, :) = [0 0 0 x y 1 -x*yp -y*yp -yp];
    end
    [~, ~, V] = svd(A);
    H_manual = reshape(V(:, end), [3, 3])';

    % Homographie MATLAB
    tform = fitgeotrans(p1, p_sorted, 'projective');
    H_matlab = tform.T;

    % Calcul des points projetés
    p1_h = [p1, ones(size(p1, 1), 1)];
    p_manual_h = (H_manual * p1_h')';
    p_matlab_h = (H_matlab * p1_h')';
    p_manual = p_manual_h(:, 1:2) ./ p_manual_h(:, 3);
    p_matlab = p_matlab_h(:, 1:2) ./ p_matlab_h(:, 3);

    % Calcul des erreurs
    error_manual = sqrt(sum((p_sorted - p_manual).^2, 2));
    error_matlab = sqrt(sum((p_sorted - p_matlab).^2, 2));
    errors_manual = [errors_manual; error_manual];
    errors_matlab = [errors_matlab; error_matlab*10e-4];
end

% Affichage des résultats
fprintf('\nErreurs moyennes :\n');
fprintf('Homographie manuelle : %.4f pixels\n', mean(errors_manual));
fprintf('Homographie MATLAB : %.4f pixels\n', mean(errors_matlab));

% Courbes d'erreur
figure;
plot(errors_manual, '-o', 'DisplayName', 'Manuelle');
hold on;
plot(errors_matlab, '-x', 'DisplayName', 'MATLAB');
legend('show');
xlabel('Index des points');
ylabel('Erreur (pixels)');
title('Comparaison des erreurs d\homographie');