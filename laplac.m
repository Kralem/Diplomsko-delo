function laplac(subj, rec)
%funkcija za parametre prebere ime direktorija ter številko prvega posnetka
%opomba: ne želimo brati prvih dveh posnetkov, saj te nimajo aktivnosti, ki
%jih iščemo
  db="eegmmidb";
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
    izpis=strcat("Racunam laplaca za subjekta ",subject);
    disp(izpis);
    for n=1:length(sig(1,:)) %laplacova formula, pdf št. 8
        sig_tmp(1,n) = sig(13,n) - 1/4 * (sig(11,n) + sig(36,n) + sig(42,n) + sig(53,n));
        sig_tmp(2,n) = sig(9,n) - 1/4 * (sig(11,n) + sig(32,n) + sig(41,n) + sig(49,n));
    end
    sig=sig_tmp;
    izpis=strcat("Berem intervale za subjekta ",subject);
    disp(izpis);
    for j=1:size(T1, 1)  %preberemo intervale zamišljanja motoričnih aktivnosti
      t1s{end+1}=sig(:, T1(j,1):T1(j,2));
    end
    for j=1:size(T2,1)
      t2s{end+1}=sig(:, T2(j,1):T2(j,2));
    end
  end
  f = [0 8 8 13 13 fs/2]/(fs/2); %definicija filtra
  a = [0 0 1 1 0 0 ];
  n = 35; 
  b = firls(n, f, a);
  
  t1s(1)=[];
  t2s(1)=[];
 
  lvt1=[];
  lvt2=[];
 
  for i=1:size(t1s,2)
      tmp = cell2mat(t1s(i));
      tmp = [tmp(1,:).',tmp(size(tmp,1),:).'].';
      tmp = filter(b,1,tmp);
      lvt1(i,1) = log(var(tmp(1,:)));
      lvt1(i,2) = log(var(tmp(2,:)));      
  end
  
  for i=1:size(t2s,2)
      tmp = cell2mat(t2s(i));
      tmp = [tmp(1,:).',tmp(size(tmp,1),:).'].';
      tmp = filter(b,1,tmp);
      lvt2(i,1) = log(var(tmp(1,:)));
      lvt2(i,2) = log(var(tmp(2,:)));      
  end
  dimenzija=size(sig);
  dimenzija %izpis dimenzije, mora biti približno 2x20000
  
  scatter(lvt1(:,1), lvt1(:,2)); %diagram raztrosa
  hold on
  scatter(lvt2(:,1), lvt2(:,2)); 
  
  featVFile = strcat(subject,'featureVectorsL.txt');
  classFile = strcat(subject,'referenceClassL.txt');
  
  fvf = fopen(featVFile, "wt");
  rcf = fopen(classFile, "wt");
  
  for i=1:size(lvt1,1)
      fprintf(fvf, "%.8f %.8f\n", lvt1(i,1), lvt1(i,2));
      fprintf(rcf, "T1\n");
  end
  
  for i=1:size(lvt2,1)
      fprintf(fvf, "%.8f %.8f\n", lvt2(i,1), lvt2(i,2));
      fprintf(rcf, "T2\n");
  end
  fclose(fvf);
  fclose(rcf);
  
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