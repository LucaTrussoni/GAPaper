% Applicazione 1 dell'articolo. Algoritmo genetico per selezionare un
% campione normale di 100 elemento da una distribuzione di 10000 per
% ottenere media e varianza ottimi
rng (6,'twister') % 1,6,18,22,33 good
X_Big=randn(10000,1);
Popolazione=ceil(rand(100,500)*10000);
% make sure all are individuals
for i=1:500
    while length(unique(Popolazione(:,i)))<100
        Popolazione(:,i)=ceil(rand(100,1)*10000);
    end
end
% 500 individui da 100 geni
MediePopolazione=mean(X_Big(Popolazione));
DevStPopolazione=std(X_Big(Popolazione));
SkewPopolazione=skewness(X_Big(Popolazione));
KurtPopolazione=kurtosis(X_Big(Popolazione));
w1=1500;
w2=10000;
w3=500;
w4=250;
Fitness=w1*MediePopolazione.^2+w2*(DevStPopolazione-1).^2+w3*SkewPopolazione.^2+w4*(KurtPopolazione-3).^2;
pdeath=0.1;
pfight=0.05;
pmut=0.5;
pcross=0.35;
genimut=5;
genicross=5;
Res=zeros(20001,1);
for i=1:10000 % generazioni
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
    elseif caso<pdeath+pfight+pmut % step (d)
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
        % crossover, il candidato � l'elemento i target; step (e)
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
    MediePopolazione=mean(X_Big(Popolazione));
    DevStPopolazione=std(X_Big(Popolazione));
    SkewPopolazione=skewness(X_Big(Popolazione));
    KurtPopolazione=kurtosis(X_Big(Popolazione));
    Fitness=w1*MediePopolazione.^2+w2*(DevStPopolazione-1).^2+w3*SkewPopolazione.^2+w4*(KurtPopolazione-3).^2;
    if (Fitness(1)~=Res(i))
        disp(Fitness(1));
    end
    if (mod(i,100)==0)
        disp(i)
    end
    Res(i+1)=Fitness(1);
end

subplot(1,2,1)
normplot(X_Big(Popolazione(:,1)))
subplot(1,2,2)
histfit(X_Big(Popolazione(:,1)),15)
