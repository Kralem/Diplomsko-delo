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
    if i~=89 && i ~= 88 && i ~= 92 && i~=100 && i~=104 %problematični subjekti
        laplac(out, 4)
        movefile(strcat(out,'featureVectorsL.txt'),'rezultatiL'); %premaknemo .txt datoteki
        movefile(strcat(out,'referenceClassL.txt'),'rezultatiL'); %na levi ime .txt datoteke, na desni ime direktorija
    end
end