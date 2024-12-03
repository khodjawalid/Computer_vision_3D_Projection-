
vectmatr = []
vectmatb = []
vectmatg = []
vectmaty = []

for i = 1:435
I = imread(sprintf('myVideoFrames2/000%03d.jpg',i));
[BWb,Ib] = createMaskb(I);
[BWr,Ir] = createMaskr(I);
[BWg,Ig] = createMaskg(I);
[BWy,Iy] = createMasky(I);


Db=BWb;
[Ly,numy]=bwlabel(Db,4); % etiquetage région 4 connexité , Ly contient les regions  
Hy=hist(Ly(:),[0:numy]); 
    for j=2:length(Hy)
        if Hy(j)<60 % on a choisi le seuil 32 pour determiner les petites regions
            Ly(Ly==j-1)=0; %suppression des petite regions, Ly==j-1 récuprère l'etiquete j-1
        else
            [yy,xy]=find(Ly==j-1);
            vectmatb=[vectmatb; [mean(xy),mean(yy)]]; %calcul du centre de gravité des région et ajout au vecteur   
        end 
    end
    

Db=BWr;
[Ly,numy]=bwlabel(Db,4); % etiquetage région 4 connexité , Ly contient les regions  
Hy=hist(Ly(:),[0:numy]); 
    for j=2:length(Hy)
        if Hy(j)<50 % on a choisi le seuil 32 pour determiner les petites regions
            Ly(Ly==j-1)=0; %suppression des petite regions, Ly==j-1 récuprère l'etiquete j-1
        else
            [yy,xy]=find(Ly==j-1);
            vectmatr=[vectmatr; [mean(xy),mean(yy)]]; %calcul du centre de gravité des région et ajout au vecteur   
        end 
    end
Db=BWg;
[Ly,numy]=bwlabel(Db,4); % etiquetage région 4 connexité , Ly contient les regions  
Hy=hist(Ly(:),[0:numy]); 
    for j=2:length(Hy)
        if Hy(j)<70 % on a choisi le seuil 32 pour determiner les petites regions
            Ly(Ly==j-1)=0; %suppression des petite regions, Ly==j-1 récuprère l'etiquete j-1
        else
            [yy,xy]=find(Ly==j-1);
            vectmatg=[vectmatg; [mean(xy),mean(yy)]]; %calcul du centre de gravité des région et ajout au vecteur   
        end 
    end
Db=BWy;
[Ly,numy]=bwlabel(Db,4); % etiquetage région 4 connexité , Ly contient les regions  
Hy=hist(Ly(:),[0:numy]); 
    for j=2:length(Hy)
        if Hy(j)<50 % on a choisi le seuil 32 pour determiner les petites regions
            Ly(Ly==j-1)=0; %suppression des petite regions, Ly==j-1 récuprère l'etiquete j-1
        else
            [yy,xy]=find(Ly==j-1);
            vectmaty=[vectmaty; [mean(xy),mean(yy)]]; %calcul du centre de gravité des région et ajout au vecteur   
        end 
    end

end  
vectmatr = vectmatr(1:435,:)
vectmatb = vectmatb(1:435,:)
vectmatg = vectmatg(1:435,:)
vectmaty = vectmaty(1:435,:)

% tracés des trajectoire

figure(11),hold on,colormap,title('trajectoire calculée avec le masque matlab du centre de gravité jaune');
 plot(vectmaty(:,1),vectmaty(:,2), '+r');
 
 figure(12),hold on,colormap,title('trajectoire calculée avec le masque matlab du centre de gravité rouge');
 plot(vectmatr(:,1),vectmatr(:,2), '*r');
 
  figure(13),hold on,colormap,title('trajectoire calculée avec le masque matlab du centre de gravité blue');
 plot(vectmatb(:,1),vectmatb(:,2), '+b');
  
  figure(14),hold on,colormap,title('trajectoire calculée avec le masque matlab du centre de gravité vert');
 plot(vectmatg(:,1),vectmatg(:,2), '+g');
 

errb = vectb - vectmatb;
errg = vectg - vectmatg;
errR = vectr - vectmatr;
erry = vecty - vectmaty;

for i = 1:435
    errb1(i) = norm(errb(i));
errg1(i) = norm(errg(i,:));
errR1(i) = norm(errR(i,:));
erry1(i) = norm(erry(i,:));
end

% tracés des erreurs

figure(),hold on,colormap,title('erreur en pixel sur le centre de gravité jaune');
 bar([1:150],erry1(1,1:150), 'y');
 
 figure(),hold on,colormap,title('erreur en pixel sur le centre de gravité rouge');
 bar([1:150],errR1(1,1:150), 'r');
 
  figure(),hold on,colormap,title('erreur en pixel sur le centre de gravité blue');
 bar([1:150],errb1(1,1:150), 'b');
  
  figure(),hold on,colormap,title('erreur en pixel sur le centre de gravité vert');
 bar([1:150],errg1(1,1:150), 'g');

