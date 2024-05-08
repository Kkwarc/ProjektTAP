function Sz=DMCstepmatricesZ(Tp,timefinal)
% konstrukcja macierzowej odpowiedzi skokowej zakłóceniowej dla DMC,
% zakładając timefinal= Tp*Dz, Dz -- liczba elementów odp. skokowej
%
% przykład obiektu: transmitancje ciągłe zakłóceniowe (WB):
s=tf('s');
WBtfZ1=3.8*exp(-8*s)/(14.9*s+1);
WBtfZ2=4.9*exp(-3*s)/(13.2*s+1);
WBtfZ=[WBtfZ1;WBtfZ2];
WBdtfZ=c2d(WBtfZ,Tp);
% Dyskretna odpowiedź skokowa wielowymiarowa wg konwencji Matlaba:
Ydstepz=step(WBdtfZ,timefinal);%macierz wymiaru (timefinal/Tp)+1 x ny x nz na odcinku
% czasu od t=0 do t=timefinal z krokiem Tp
% Dyskretna odpowiedź skokowa macierzowa Sz o wymiarze ny x nz x timefinal/Tp; pierwsza
% macierz dla czasu dyskr. k=1, ostatnia dla k=Dz=timefinal/Tp (Dz macierzy):
[nt,ny,nz]=size(Ydstepz); % nt=D+11 (Ydstepz zawiera wartość wyjść w chwili 0, Sz nie)
Sz=zeros(ny,nz,nt-1);
for i=1:ny
for j=1:nz
Sz(i,j,:)=Ydstepz(2:nt,i,j); %
end
end
