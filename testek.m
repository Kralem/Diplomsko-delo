s='S00';
for i=1:3
    if i >= 10
        s='S0';
    elseif i >= 100
        s='S';
    end
    
    out=strcat(s,num2str(i));
    disp(out);
    laplac(out, 4)
end