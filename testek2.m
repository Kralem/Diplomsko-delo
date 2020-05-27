s='S00';
for i=88:109
    if i >= 10
        s='S0';
    end
    if i >= 100 && i >= 10
        s='S';
    end
    
    disp(num2str(i));
    out=strcat(s,num2str(i));
    disp(out);
    if i~=89 && i ~= 88 && i ~= 92 && i~=100 && i~=104 %S089 je problematičen subjekt, zato ga izpustimo
        izracunZnacilkPopravljeno(out, 4)
        movefile(strcat(out,'featureVectors.txt'),'rezultatiC'); %premaknemo .txt datoteki
        movefile(strcat(out,'referenceClass.txt'),'rezultatiC'); %na levi ime .txt datoteke, na desni ime direktorija
    end    
end