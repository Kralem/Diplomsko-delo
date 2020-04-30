[sig, fs, tm] = rdsamp('/eegmmidb/S001/S001R08.edf', 1:64); % ker ne rabimo zadnjega kanala

[T0, T1, T2] = getIntervals('/eegmmidb/S001/S001R08.edf','event', 160);
sig=sig'; % matriko obrnemo zaradi lazjega racunanja

% v c1dat si zapomnimo vzorce T1 iz zacetka posnetka (ucenje)
% v c2dat si zapomnimo vzorce T2 iz zacetka posnetka (ucenje)
% v c3dat si zapomnimo vzorce T1 iz konca posneta (operiranje)
% v c4dat si zapomnimo vzorce T2 iz konca posnetka (operiranje)
c1dat = sig(:, T1(1,1):T1(1,2)); 
c2dat = sig(:, T2(1,1):T2(1,2));
c3dat = sig(:, T1(size(T1,1), 1):T1(size(T1,1), 2)); 
c4dat = sig(:, T2(size(T2,1), 1):T2(size(T2,1), 2));

% dodamo se tocke za T1 in T2!!!
c3dat2 = sig(:, T1(size(T1,1)-1, 1):T1(size(T1,1)-1, 2));% 16609 + 4.1*fs - T1
c4dat2 = sig(:, T2(size(T2,1)-1, 1):T2(size(T2,1)-1, 2));% 13953 + 4.1*fs - T2
c3dat3 = sig(:, T1(size(T1,1)-2, 1):T1(size(T1,1)-2, 2));% 15281 + 4.1*fs - T1
c4dat3 = sig(:, T2(size(T2,1)-2, 1):T2(size(T2,1)-2, 2));% 12625 + 4.1*fs - T2

% izrisemo signale
plot(c1dat)
% sedaj pa se pravilen izris
figure
plot(c1dat');
% odstejemo se srednjo vrednost, kar ni povsem nujno, in to izrisemo
% [c1dat mean] = remmean(c1dat);
% figure
% plot (c1dat');
% odstejemo se srednjo vrednost za ostale signale, lahko uporabljamo
% odstete vrednost ali pa neodstete
% [c2dat mean] = remmean(c2dat);
% [c3dat mean] = remmean(c3dat);
% [c4dat mean] = remmean(c4dat);

% [c3dat2 mean] = remmean(c3dat2);
% [c4dat2 mean] = remmean(c4dat2);
% [c3dat3 mean] = remmean(c3dat3);
% [c4dat3 mean] = remmean(c4dat3);

% sedaj pa lahko pridobimo matriko filtrov za maksimiziranje W iz obeh primerov 
[W] = f_CSP(c1dat, c2dat);
% sedaj izpeljemo protor komponent za oba primera, stiskanje leve oziroma
% desne roke
S3 = W*c3dat;
S4 = W*c4dat;
% dodamo se primere za stiskanje leve in desne roke, ki smo jih pridobili
% na podlagi oznak v datoteki z oznakami
S32 = W*c3dat2;
S42 = W*c4dat2;
S33 = W*c3dat3;
S43 = W*c4dat3;

% sedaj konstruiramo vektorje za c3dat ter c4dat, pri cemer uporabljamo
% samo prvi in zadnji signal v prostoru komponent, tistega z min-max
% varianco za dve stanji
ff=S3(1,:).';
fl=S3(size(S3,1), :).';
S3A=transpose([ff fl]);
S4A=[S4(1,:).', S4(size(S4,1),:).' ].';
figure
plot (S3A.');
figure
plot (S4A.');
figure
plot(S3A(1,:), S3A(2,:), 'o');
hold on;
plot (S4A(1,:), S4A(2,:),'x');
S32A = [S32(1,:).', S32(size(S32,1),:).' ].';
S42A = [S42(1,:).', S42(size(S42,1),:).' ].'; 
S33A = [S33(1,:).', S33(size(S33,1),:).' ].';
S43A = [S43(1,:).', S43(size(S43,1),:).' ].';

for i=1:2
    lv3A(i,:)  = log(var(S3A(i,:)));
    lv4A(i,:)  = log(var(S4A(i,:)));
    lv32A(i,:) = log(var(S32A(i,:)));
    lv42A(i,:) = log(var(S42A(i,:)));
    lv33A(i,:) = log(var(S33A(i,:)));
    lv43A(i,:) = log(var(S43A(i,:)));
end


function [newVectors, meanValue] = remmean(vectors);
%
% Remove the mean from vectors
%
% Removes the mean of row vectors.
% Returns the new vectors and the mean.
%
% @(#)$Id: remmean.m,v 1.2 2003/04/05 14:23:58 jarmo Exp $

newVectors = zeros (size (vectors));
meanValue = mean (vectors')';
newVectors = vectors - meanValue * ones (1,size (vectors, 2));

end
