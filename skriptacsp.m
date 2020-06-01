s='S00';
for i=1:109
    if i >= 10
        s='S0';
    end
    if i >= 100 && i >= 10
        s='S';
    end
    out=strcat(s,num2str(i));
    if i~=89 && i ~= 88 && i ~= 92 && i~=100 && i~=104 %S089 je problematičen subjekt, zato ga izpustimo
        disp(out);
        vec=strcat(out,'featureVectors.txt');
        ref=strcat(out,'referenceClass.txt');
        doClassification(vec,ref, {1,1}, 10, 50, 0);
    end    
end