function S=DMCstepmatrices(Tp,timefinal)
    % konstrukcja macierzowej odpowiedzi skokowej dla DMC, zakładając timefinal=Tp*D
    % Tp - okres próbkowania, D - horyzont dynamiki obiektu.
    % Przykład obiektu: transmitancje ciągłe obiektu 2x2 (kolumna WB):
    s=tf('s');
    WBtf11=0.48/(145.1*s+1);
    WBtf12=0.48*exp(-90*s)/(145.1*s+1);
    % WBtf21=0;
    % WBtf22=(1254*s-8.64)/(145.1*s+1);
    WBtf31=(97.2*s+0.6698)*exp(-40*s)/(3510*s^2+169.3*s+1);
    WBtf32=(47.92*s+0.3302)*exp(-130*s)/(3510*s^2+169.3*s+1);
    % WBtf41=0.1667*exp(-40*s)/(24.19*s+1);
    % WBtf42=(129300*s^2+5345*s+30.69)*exp(-40*s)/(3510*s^2+169.3*s+1);
    WBtf=[WBtf11 WBtf12; WBtf31 WBtf32];
    % WBtf=[WBtf11 WBtf12; WBtf21 WBtf22; WBtf31 WBtf32; WBtf41 WBtf42];
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
