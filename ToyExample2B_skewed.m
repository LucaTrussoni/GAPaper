% Applicazione 1 dell'articolo. Algoritmo genetico per selezionare un
% campione normale di 100 elemento da una distribuzione di 10000 per
% ottenere media e varianza ottimi
rng (1,'twister') % 2,13,32,51 good
Nper=12;
sigma=0.15;
r=0.04;
dt=1/Nper;
X_Big_Processes=r*dt+randn(10000,Nper)*sigma*(dt^0.5)-0.5*(sigma^2)*dt;
% dodici rendimenti mensili, tasso annuale 0.04, volatilit� 0.15
X_Big_Rendimenti=sum(X_Big_Processes,2);
X_Big_CashAccount=exp(X_Big_Rendimenti);
Popolazione=ceil(rand(100,500)*10000);
% make sure all are individuals
for i=1:500
    while length(unique(Popolazione(:,i)))<100
        Popolazione(:,i)=ceil(rand(100,1)*10000);
    end
end
% 1000 individui da 100 geni
MediePopolazione=mean(X_Big_CashAccount(Popolazione));
DevStPopolazione=std(X_Big_Rendimenti(Popolazione));
SkewPopolazione=skewness(X_Big_Rendimenti(Popolazione)/sigma);
KurtPopolazione=kurtosis(X_Big_Rendimenti(Popolazione)/sigma)-3;
w1=125;
w2=30;
w3=5;
w4=1;
Fitness=w1*(MediePopolazione-exp(r)).^2+w2*(DevStPopolazione-sigma).^2+w3*SkewPopolazione.^2+w4*KurtPopolazione.^2;
pdeath=0.1;
pfight=0.1;
pmut=0.4;
pcross=0.3;
Res=zeros(10001,1);
genimut=5;
genicross=5;
for i=1:15000 % generazioni
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
    elseif caso<pdeath+pfight % competizione con avversario casuale step (c)
        target=ceil(rand()*500);
        winner=icaso;
        if Fitness(icaso)<Fitness(target)
            winner=target;
        end
        winnerG=Popolazione(:,winner);
        Popolazione(:,icaso)=winnerG;
        Popolazione(:,target)=winnerG; % il genoma dei contendenti � sostituito da quello del vincitore
    elseif caso<pdeath+pfight+pmut % mutazione, step (d)
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
    SkewPopolazione=skewness(X_Big_Rendimenti(Popolazione)/sigma);
    KurtPopolazione=kurtosis(X_Big_Rendimenti(Popolazione)/sigma)-3;
    Fitness=w1*(MediePopolazione-exp(r)).^2+w2*(DevStPopolazione-sigma).^2+w3*SkewPopolazione.^2+w4*KurtPopolazione.^2;
    %Fitness=w1*(MediePopolazione-exp(r)).^2+w2*(DevStPopolazione-sigma).^2;
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
plot(cumsum([zeros(100,1),X_Big_Processes(Popolazione(:,1),:)],2)')
disp(mean(X_Big_CashAccount(Popolazione(:,1))))
disp(std(log(X_Big_CashAccount(Popolazione(:,1)))))
disp(skewness(log(X_Big_CashAccount(Popolazione(:,1)))))
disp(kurtosis(log(X_Big_CashAccount(Popolazione(:,1)))))
