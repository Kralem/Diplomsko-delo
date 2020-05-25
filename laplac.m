function laplac(subj, rec)
%funkcija za parametre prebere ime direktorija ter številko prvega posnetka
%opomba: ne želimo brati prvih dveh posnetkov, saj te nimajo aktivnosti, ki
%jih iščemo
  db="baza";
  subject = string(subj);
  record=[];
  for i=0:2
      record = [record string(num2str(rec+(i*4), '%02d'))];
  end
  t1s = {};
  t2s = {};
  
  sig_tmp = []; %definiramo sig_tmp, ki ga uporabimo pri laplacovi maski
  
  for i=1:size(record,2)
    if (strcmp(db,"")==0) 
      recName=strcat("/",db,"/",subject,"/",subject,"R",record(:,i),".edf");
    else
      recName=strcat(subject,'R',record(i),'.edf');
    end
    recName=convertStringsToChars(recName);
    disp(recName);
    [sig, fs, tm] = rdsamp(recName, 1:64); %branje posnetka
    [T0, T1, T2] = getIntervals(recName,'event', fs, size(sig,1)); %pridobitev intervalov
    sig=sig';
    izpis=strcat("Racunam laplaca za subjekta ",recName);
    disp(izpis);
    for n=1:length(sig(1,:)) %laplacova formula, pdf št. 8
        sig_tmp(1,n) = sig(13,n) - 1/4 * (sig(11,n) + sig(36,n) + sig(42,n) + sig(53,n));
        sig_tmp(2,n) = sig(9,n) - 1/4 * (sig(11,n) + sig(32,n) + sig(41,n) + sig(49,n));
    end
    
    sig=sig_tmp;
    izpis=strcat("Berem intervale za subjekta ",recName);
    disp(izpis);
    for j=1:size(T1, 1)  %preberemo intervale zamišljanja motoričnih aktivnosti
      t1s{end+1}=sig(:, T1(j,1):T1(j,2));
    end
    for j=1:size(T2,1)
      t2s{end+1}=sig(:, T2(j,1):T2(j,2));
    end
  end
  
  if size(t1s,2) > 0 && size(t2s,2) > 0 %varovalka, če subjekt nima intervalov zamišljanja T1 ali T2
      
      
      %figure;
      %plot(cell2mat(t1s(1,:)));
      %figure;
      %plot(cell2mat(t2s(:,1)));
      f = [0  7  8  13  14  fs/2]/(fs/2);  %definicija filtra
      a = [0 0 1 1 0 0 ];
      n = 35; 
      b = firls(n, f, a);
      
      %win_index_shift = 100; % ce ne vizualizira signale zamisljanj obeh rok v isto sliko
      
      %t1s(1)=[];
      %t2s(1)=[];

      lvt1=[];
      lvt2=[];

lengths_of_intervals_in_samlpes_left = [];
lengths_of_intervals_in_seconds_left = [];
%%%%

  for i=1:size(t1s,2)
      % tmp = (W*filter(b, 1,cell2mat(t1s(i))));

% Signali se prevedejo v prostor stanj 64 -> 64,
% zanimata nas le prvi (postane PRVI signal v prostoru stanj) 
% in 64-ti (postane DRUGI signal v prostoru stanj)

        tmp = (cell2mat(t1s(i)));

% disp(length(tmp(1,:)));
% disp(length(tmp(64,:)));
% pause

% Slike signalov za intervale zamisljanja leve roke pred filtriranjem (slike od 1 naprej)

%figure(i)

%%%%     length(tmp(1,:)))       dolzina zamisljanja aktivnosti leve roke v vzorcih za ta interval
%%%%     length(tmp(1,:))*1/fs)  dolzina zamisljanja aktivnosti leve roke v sekundah za ta interval

%%%%    tmp(1,:)                 prvi signal v prostoru stanj (na polozaju 1 v matriki) pred filtriranjem
%plot(tmp(1,:))                 % Prvi signal, Slika je v samplih signala po X osi
%ylim([-1000, 1000]);
%hold on                        % naslednji signa bo v isti sliki
%plot(tmp(size(tmp,1),:))       % Drugi signal, Slika je v samplih signala po X osi
%ylim([-1000, 1000]);
%%%%    tmp(size(tmp,1),:)       drugi signal v prostoru stanj (na polozaju 64 v matriki) pred filtriranjem

lengths_of_intervals_in_samlpes_left = [lengths_of_intervals_in_samlpes_left, length(tmp(1,:))];
lengths_of_intervals_in_seconds_left = [lengths_of_intervals_in_seconds_left, length(tmp(1,:))*1/fs];

% Tu se filtrirata prvi in drugi signal v prostoru stanj
      tmp = [tmp(1,:).', tmp(size(tmp,1),:).' ].';
      tmp = filter(b, 1, tmp);

% Slike signalov za intervale zamisljanja leve roke PO filtriranju (slike od 100 naprej)

%figure(i+win_index_shift)
%%%%     length(tmp(1,:)))       dolzina zamisljanja aktivnosti leve roke v vzorcih za ta interval
%%%%     length(tmp(1,:))*1/fs)  dolzina zamisljanja aktivnosti leve roke v sekundah za ta interval

%%%%     tmp(1,:)                prvi signal v prostoru stanj PO filtriranju (na polozaju 1 v matriki)
%plot(tmp(1,:))  % Prvi signal,   Slika je v samplih signala po X osi
%ylim([-100, 100]);
%hold on
%plot(tmp(2,:))  % Drugi signal,  Slika je v samplih signala po X osi
%ylim([-100, 100]);
%%%%     tmp(2,:)                drugi signal v prostoru stanj PO filtriranju (na polozaju 2 v matriki)

%% Na tem mestu se filtrirani intervali zamisljanj leve roke (po dva signala v prostoru stanj)
%% lahko zapisejo na disk

% Izracun znacilk za intervale zamisljanja aktivnosti leve roke
      lvt1(i,1) = log(var(tmp(1,:)));
      lvt1(i,2) = log(var(tmp(2,:)));     
  end

[mrows, ncols] = size(lengths_of_intervals_in_seconds_left);
outputstr = ['%' num2str(mrows) '.4f ']; % template for the string, you put your datatype here
outputstr = repmat(outputstr, 1, ncols); % replicate it to match the number of columns
outputstr = [outputstr '\n']; % add a new line if you want 
  
file = 'intervali.txt';
fp = fopen(file, "at+");
fprintf(fp, "Subjekt: %s \n", subj);
fprintf(fp, outputstr, lengths_of_intervals_in_seconds_left);
fprintf(fp, "T1: %d \n", size(lvt1,1));
  

% Dolzine intervalov zamisljanja aktivnosti leve roke v vzorcih
lengths_of_intervals_in_samlpes_left
% Dolzine intervalov zamisljanja aktivnosti leve roke v sekundah
lengths_of_intervals_in_seconds_left

%%%%

lengths_of_intervals_in_samlpes_right = [];
lengths_of_intervals_in_seconds_right = [];

  for i=1:size(t2s,2)
      % tmp = (W*filter(b, 1,cell2mat(t2s(i))));

% Signali se prevedejo v prostor stanj 64 -> 64,
% zanimata nas le prvi (postane PRVI signal v prostoru stanj) 
% in 64-ti (postane DRUGI signal v prostoru stanj)

        tmp = (cell2mat(t2s(i)));

% disp(length(tmp(1,:)));
% disp(length(tmp(64,:)));
% pause

% Slike signalov za intervale zamisljanja desne roke pred filtriranjem (slike od 200 naprej)

%figure(i+2*win_index_shift)

%%%%    length(tmp(1,:)))        dolzina zamisljanja aktivnosti desne roke v vzorcih za ta interval
%%%%    length(tmp(1,:))*1/fs)   dolzina zamisljanja aktivnosti desne roke v sekundah za ta interval

%%%%    tmp(1,:)                 prvi signal v prostoru stanj (na polozaju 1 v matriki) pred filtriranjem
%plot(tmp(1,:))                 % Prvi signal, Slika je v samplih signala po X osi
%ylim([-1000, 1000]);
%hold on                        % naslednji signa bo v isti sliki
%plot(tmp(size(tmp,1),:))       % Drugi signal, Slika je v samplih signala po X osi
%ylim([-1000, 1000]);
%%%%    tmp(size(tmp,1),:)       drugi signal v prostoru stanj (na polozaju 64 v matriki) pred filtriranjem

lengths_of_intervals_in_samlpes_right = [lengths_of_intervals_in_samlpes_right, length(tmp(1,:))];
lengths_of_intervals_in_seconds_right = [lengths_of_intervals_in_seconds_right, length(tmp(1,:))*1/fs];

% Tu se filtrirata prvi in drugi signal v prostoru stanj
      tmp = [tmp(1,:).', tmp(size(tmp,1),:).' ].';
      tmp = filter(b, 1, tmp);

% Slike signalov za intervale zamisljanja desne roke PO filtriranju (slike od 300 dalje)

%figure(i+3*win_index_shift)

%%%%    length(tmp(1,:)))        dolzina zamisljanja aktivnosti desne roke v vzorcih za ta interval
%%%%    length(tmp(1,:))*1/fs)   dolzina zamisljanja aktivnosti desne roke v sekundah za ta interval

%%%%    tmp(1,:)                 prvi signal v prostoru stanj PO filtriranju (na polozaju 1 v matriki)
%plot(tmp(1,:))  % Prvi signal,   Slika je v samplih signala po X osi
%ylim([-100, 100]);
%hold on                        % naslednji signa bo v isti sliki
%plot(tmp(2,:))  % Drugi signal,  Slika je v samplih signala po X osi
%ylim([-100, 100]);
%%%%    tmp(2,:)                 drugi signal v prostoru stanj PO filtriranju (na polozaju 2 v matriki)

%% Na tem mestu se filtrirani intervali zamisljanj leve roke (po dva signala v prostoru stanj)
%% lahko zapisejo na disk

% Izracun znacilk za intervale zamisljanja aktivnosti desne roke roke
      lvt2(i,1) = log(var(tmp(1,:)));
      lvt2(i,2) = log(var(tmp(2,:)));     
  end

[mrows2, ncols2] = size(lengths_of_intervals_in_seconds_right);
outputstr2 = ['%' num2str(mrows2) '.4f ']; % template for the string, you put your datatype here
outputstr2 = repmat(outputstr2, 1, ncols2); % replicate it to match the number of columns
outputstr2 = [outputstr2 '\n']; % add a new line if you want   
  
fprintf(fp, outputstr2, lengths_of_intervals_in_seconds_right);
fprintf(fp, "T2: %d \n", size(lvt2,1));
fprintf(fp, "================================================================== \n");
fclose(fp);

% Dolzine intervalov zamisljanja aktivnosti desne roke v vzorcih
lengths_of_intervals_in_samlpes_right
% Dolzine intervalov zamisljanja aktivnosti desne roke v sekundah
lengths_of_intervals_in_seconds_right

%%%



%figure(4*win_index_shift)  %% Slika 500 Diagram raztrosa znacilk

      %scatter(lvt1(:,1), lvt1(:,2));
      %hold on
      %scatter(lvt2(:,1), lvt2(:,2)); 
     
      dimenzija=size(sig);
      dimenzija %izpis dimenzije, mora biti približno 2x20000


      %featVFile = strcat(subject,'featureVectorsL.txt');
      %classFile = strcat(subject,'referenceClassL.txt');
      %fvf = fopen(featVFile, "wt");
      %rcf = fopen(classFile, "wt");

      %for i=1:size(lvt1,1)
      %    fprintf(fvf, "%.8f %.8f\n", lvt1(i,1), lvt1(i,2));
      %    fprintf(rcf, "T1\n");
      %end

      %for i=1:size(lvt2,1)
      %    fprintf(fvf, "%.8f %.8f\n", lvt2(i,1), lvt2(i,2));
      %    fprintf(rcf, "T2\n");
      %end
      %fclose(fvf);
      %fclose(rcf);
  
  end
  
end
  
function [T0,T1, T2] = getIntervals(recName, anotator, fs, siglen)
% function returns arrays of interval start and stop times in samples 
% for the three different interval types; the meaning of annotations (T1,
% T2, T3) is determined according to the record number;
% input arguments:
%   - recName: string containing record name e.g 'S001R03.edf'
%   - anotator: string containing annotator name e.g. 'event'
%   - fs: sampling frequency of the record e.g. 160
%   - siglen: length of a data file in samples
% output arguments:
%   - T0: array of start and end times in samples for interval T0
%   - T1: array of start and end times in samples for interval T1
%   - T2: array of start and end times in samples for interval T2
  T0=[];
  T1=[];
  T2=[];
  [a, t, s, c, n, cmt] = rdann(recName, anotator);
  cmt=string(cmt);
  for i=1:size(a,1)
    splt=strsplit(cmt(i), " ");
    dur=str2double(splt(3))*fs;
    if (i<size(a,1)) 
        stop=a(i+1)-1;
    else
        stop=a(i)+dur;
    end
    stop = min(stop, siglen);
    if (splt(1)=="T0") 
      T0=[T0; [a(i) stop]];
    elseif (splt(1)=="T1")
      T1=[T1; [a(i) stop]];
    elseif (splt(1)=="T2")
      T2=[T2; [a(i) stop]];
    end
  end
end

function [newVectors, meanValue] = remmean(vectors)
  newVectors = zeros (size (vectors));
  meanValue = mean (vectors')';
  newVectors = vectors - meanValue * ones (1,size (vectors, 2));
end