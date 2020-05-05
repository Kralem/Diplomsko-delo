s='S00';
for i=1:5
    if i >= 10
        s='S0';
    elseif i >= 100
        s='S';
    end
    
    out=strcat(s,num2str(i));
    disp(out);
    laplac(out, 4) %kličemo funkcijo laplace za izračun značilk
    movefile(strcat(out,'featureVectorsL.txt'),'rezultatiL'); %premaknemo .txt datoteki
    movefile(strcat(out,'referenceClassL.txt'),'rezultatiL'); %na levi ime .txt datoteke, na desni ime direktorija
end