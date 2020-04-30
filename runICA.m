%izracun sinteticnih signalov
[ sig , mix ] = demosig () ;
% izris sinteticnih komponent
h = zeros(size(sig,1));
for i =1: size ( sig ,1)
  h ( i ) = subplot ( size ( sig ,1) ,1 , i ) ;
  plot ( sig (i ,:) ) ;
end
%izris sinteticnih signalov
figure
k = zeros(size(sig,1));
for i =1: size ( mix ,1)
  k ( i ) = subplot ( size ( mix ,1) ,1 , i ) ;
  plot ( mix (i ,:) ) ;
end

%branje dejanskih signalov
%bazalni nivo, odprte oci
%[sig, freq, tm] = rdsamp('S002R01.edf');
%plot (tm, sig(:,1:64));
%figure;
%bazalni nivo, zaprte oci
%[sig1, freq, tm] = rdsamp('S002R02.edf');
%plot (tm, sig1(:,1:64));
%figure

%branje samo stirih signalov, zaradi enostavnosti prikaza
%intvl=1:1600;
%[sig, freq, tm] = rdsamp('S002R01.edf',1:4);
%sig=sig';
%sig=sig(:, intvl);
%h=[];
%for i=1:size(sig,1)
%  h(i)=subplot(size(sig,1),1,i);
%  plot(tm(intvl), sig(i,:));
%end



[ icasig , A , W ] = fastica ( mix ) ;
% za prave signale potem vzamemo drugo matriko signalov
%[ icasig , A , W ] = fastica ( sig ) ;


W
A
inv ( W )
A - inv ( W )

figure
l = zeros(size(sig,1));
for i =1: size ( icasig ,1)
  l ( i ) = subplot ( size ( icasig ,1) ,1 , i ) ;
  plot ( icasig (i ,:) ) ;
end

W1 = inv ( W ) ; 
Wap = W1 (: , 2:4) ;
A (: ,2:4) - Wap 
sigs = icasig (2:4 ,:) ;
Xap = Wap * sigs ;

figure
m = zeros(size(sig,1));
for i =1: size ( Xap ,1)
  m ( i ) = subplot ( size ( Xap ,1) ,1 , i ) ;
  plot ( Xap (i ,:) ) ;
end

