function S=DMCstepmatrices(Tp,timefinal)
    % konstrukcja macierzowej odpowiedzi skokowej dla DMC, zakładając timefinal=Tp*D
    % Tp - okres próbkowania, D - horyzont dynamiki obiektu.
    % Przykład obiektu: transmitancje ciągłe obiektu 2x2 (kolumna WB):
    s=tf('s');
    WBtf11=12.8*exp(-s)/(16.7*s+1);
    WBtf12=-18.9*exp(-3*s)/(21*s+1);
    WBtf21=6.6*exp(-7*s)/(10.9*s+1);
    WBtf22=-19.4*exp(-3*s)/(14.4*s+1);
    WBtf=[WBtf11 WBtf12;WBtf21 WBtf22];
    % dyskretyzacja z okresem próbkowania Tp:
    WBdtf=c2d(WBtf,Tp);
    % Dyskretna odpowiedź skokowa wielowymiarowa wg konwencji Matlaba:
    Y=step(WBdtf,timefinal); % macierz o wymiarze (timefinal/Tp+1) x ny x nu na odcinku
    % czasu od t=0 do t=timefinal z krokiem Tp
    % Dyskretna odpowiedź skokowa macierzowa S o wymiarze ny x nu x timefinal/Tp, pierwsza
    % macierz dla czasdyskr.k=1, ostatnia dla czasu dyskr. k=D=timefinal/Tp (D macierzy):
    [nt,ny,nu]=size(Y); % nt=D+1 (Y zawiera też wartość wyjść w chwili 0, macierz S nie)
    S=zeros(ny,nu,nt-1);
    for i=1:ny
        for j=1:nu
            S(i,j,:)=Y(2:nt,i,j);
        end
    end
end
