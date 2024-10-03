% Applicazione 1 dell'articolo. Algoritmo genetico per selezionare un
% campione normale di 100 elemento da una distribuzione di 10000 per
% ottenere media e varianza ottimi
rng(56,'twister') % 56 good
Nper=12;
sigma=0.25;
a=3.5;
TimeLength=4;
dt=TimeLength/Nper;
X_Big_Random=randn(10000,Nper);
X_Big=zeros(10000,1);
cfA=exp(-a*dt);
cfB=sigma*((1-exp(-2*a*dt))/(2*a))^0.5;
for i=1:Nper
    X_Big=[X_Big,X_Big(:,i)*cfA+X_Big_Random(:,i)*cfB];
end
% dodici rendimenti mensili, tasso annuale 0.04, volatilit� 0.15
X_Big_CashAccount=exp(sum(X_Big,2)*dt);
X_Big_Rendimenti=log(X_Big_CashAccount);
sigma_CA=sigma/a*(TimeLength-2*(1-exp(-a*TimeLength))/a+(1-exp(-2*a*TimeLength))/(2*a))^0.5;
mean_CA=exp(0.5*(sigma_CA^2));
Popolazione=ceil(rand(100,500)*10000);
% make sure all are individuals
for i=1:500
    while length(unique(Popolazione(:,i)))<100
        Popolazione(:,i)=ceil(rand(100,1)*10000);
    end
end
% 500 individui da 100 geni
MediePopolazione=mean(X_Big_CashAccount(Popolazione));
DevStPopolazione=std(X_Big_Rendimenti(Popolazione),1);
SkewPopolazione=skewness(X_Big_Rendimenti(Popolazione)/sigma_CA);
KurtPopolazione=kurtosis(X_Big_Rendimenti(Popolazione)/sigma_CA)-3;
w1=50;
w2=100;
w3=10;
w4=1;
Fitness=w1*(MediePopolazione-mean_CA).^2+w2*(DevStPopolazione-sigma_CA).^2+w3*SkewPopolazione.^2+w4*KurtPopolazione.^2;
%w1*(MediePopolazione-mean_CA).^2+w2*(DevStPopolazione-sigma_CA).^2;
pdeath=0.15;
pfight=0.1;
pmut=0.35;
pcross=0.35;
genimut=5;
genicross=5;
Res=zeros(20001,1);
for i=1:20000 % generazioni
    % Salva il migliore mettendolo in prima posizione
    [a,b]=sort(Fitness);
    aux=Popolazione(:,b(1));
    auxF=Fitness(b(1));
    Popolazione(:,b(1))=Popolazione(:,1);
    Fitness(b(1))=Fitness(1);
    Popolazione(:,1)=aux;
    Fitness(1)=auxF;
    % estrae il valore caso e il valore icaso, che rappresentano cosa
    % avvrr� (death/fight/mutation/crossover) a chi (icaso). Corrisponde
    % a moltiplicare per 1/N=1/500 la probabilit� degli eventi
    caso=rand();
    icaso=ceil(rand()*499)+1; % step (a)
    if caso<pdeath % death, step (b)
        Popolazione(:,icaso)=Popolazione(:,1); % il genoma muore e viene sostituito dal migliore
    elseif caso<pdeath+pfight % competizione con avversario casuale step(c)
        target=ceil(rand()*500);
        winner=icaso;
        if Fitness(icaso)<Fitness(target)
            winner=target;
        end
        winnerG=Popolazione(:,winner);
        Popolazione(:,icaso)=winnerG;
        Popolazione(:,target)=winnerG; % il genoma dei contendenti � sostituito da quello del vincitore
    elseif caso<pdeath+pfight+pmut % Mutazione, step (d)
        for j=1:genimut
            goodmut=false;
            while not(goodmut)
                target=ceil(rand()*10000); % valore del gene mutato
                goodmut=not(ismember(target,Popolazione(:,icaso))); % check ammissibilit� mutazione
            end
            % Il target � un nuovo elemento
            Popolazione(ceil(rand()*100),icaso)=target; % la mutazione � eseguita quando ammissibile
        end
    else
        % crossover, il candidato � l'elemento i target step (e)
        goodtarget=false;
        while not(goodtarget)
            targeti=ceil(rand()*499)+1;
            goodtarget=(targeti~=icaso);
        end % ora siamo sicuri che itarget e icaso sono elementi diversi, e nessuno � il migliore
        for j=1:genicross
            goodcross=false;
            % cerchiamo genicross geni da scambiare. Due geni sono buoni da
            % scambiare se:
            % 1-sono uguali (e non c'� nulla da scambiare) oppure
            % 2-il gene dell'elemento A non compare gi� nel genoma di B e
            % viceversa
            while not(goodcross)
                targetg=ceil(rand()*100);
                if Popolazione(targetg,targeti)==Popolazione(targetg,icaso)
                    goodcross=true;
                else
                    if not(ismember(Popolazione(targetg,targeti),Popolazione(:,icaso))) && ...
                            not(ismember(Popolazione(targetg,icaso),Popolazione(:,targeti)))
                        aux=Popolazione(targetg,targeti);
                        Popolazione(targetg,targeti)=Popolazione(targetg,icaso);
                        Popolazione(targetg,icaso)=aux;
                        goodcross=true;
                    end
                end
            end
        end
    end
    MediePopolazione=mean(X_Big_CashAccount(Popolazione));
    DevStPopolazione=std(X_Big_Rendimenti(Popolazione));
    SkewPopolazione=skewness(X_Big_Rendimenti(Popolazione)/sigma_CA);
    KurtPopolazione=kurtosis(X_Big_Rendimenti(Popolazione)/sigma_CA)-3;
    Fitness=w1*(MediePopolazione-mean_CA).^2+w2*(DevStPopolazione-sigma_CA).^2+w3*SkewPopolazione.^2+w4*KurtPopolazione.^2;
    if (Fitness(1)~=Res(i))
        disp(Fitness(1));
    end
    if (mod(i,100)==0)
        disp(i)
    end
    Res(i+1)=Fitness(1);
end

subplot(1,2,1)
normplot(log(X_Big_CashAccount(Popolazione(:,1))))
subplot(1,2,2)
histfit(log(X_Big_CashAccount(Popolazione(:,1))),15)
figure()
plot(X_Big(Popolazione(:,1),:)')
disp(mean(X_Big_CashAccount(Popolazione(:,1))))
disp(std(log(X_Big_CashAccount(Popolazione(:,1)))))
disp(skewness(X_Big_Rendimenti(Popolazione(:,1))/sigma_CA))
disp(kurtosis(X_Big_Rendimenti(Popolazione(:,1))/sigma_CA))
