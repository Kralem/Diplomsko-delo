s='S00';
for i=1:109
    if i >= 10
        s='S0';
    end
    if i >= 100 && i >= 10
        s='S';
    end
    
    disp(num2str(i));
    out=strcat(s,num2str(i));
    disp(out);
    laplac(out, 4) %kličemo funkcijo laplace za izračun značilk
    if i~= 89 %S089 je problematičen subjekt, zato ga izpustimo
        movefile(strcat(out,'featureVectorsL.txt'),'rezultatiL'); %premaknemo .txt datoteki
        movefile(strcat(out,'referenceClassL.txt'),'rezultatiL'); %na levi ime .txt datoteke, na desni ime direktorija
    end
    end