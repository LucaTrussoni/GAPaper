% Example 3+: same as example 3 (GA2b.m) but we introduce
% 3rd and 4th moment in the fitness. We will need more epochs
rng(56,'twister')
Nper=12;
sigma=0.25;
a=3.5;
TimeLength=4;
dt=TimeLength/Nper;
% Creation of the pool, the normal random variable used for simulations
X_Big_Random=randn(10000,Nper);
X_Big=zeros(10000,1);
% Coefficient useful for simulation. Basically if we want to simulate
% dx=-a x dt + sigma dW in discrete time we can write
% x(t+dt)=cfA x(t) + cfB * N(0,1)
cfA=exp(-a*dt);
cfB=sigma*((1-exp(-2*a*dt))/(2*a))^0.5;
% In-place cretion of simulated paths
for i=1:Nper
    X_Big=[X_Big,X_Big(:,i)*cfA+X_Big_Random(:,i)*cfB];
end
% Creation of pools of terminal results of the simulations
X_Big_CashAccount=exp(sum(X_Big,2)*dt);
X_Big_Yields=log(X_Big_CashAccount);
% Expected values for cash account and volatility of terminal results of
% simulations (see Brigo-Mercurio book)
sigma_CA=sigma/a*(TimeLength-2*(1-exp(-a*TimeLength))/a+(1-exp(-2*a*TimeLength))/(2*a))^0.5;
mean_CA=exp(0.5*(sigma_CA^2));
% Create a population of 500 individual. Each individual has a genoma of
% 100 genes, any of them points to a path in the pool
Population=ceil(rand(100,500)*10000);
% We make sure all individuals are well defined, so no path is selected
% twice
for i=1:500
    while length(unique(Population(:,i)))<100
        Population(:,i)=ceil(rand(100,1)*10000);
    end
end
% Initialize fitness: our goal is to find an individual whose statistics
% are close to theoretical values of mean and variance
MeanPopulation=mean(X_Big_CashAccount(Population));
StDevPopulation=std(X_Big_Yields(Population),1);
SkewPopulation=skewness(X_Big_Yields(Population)/sigma_CA,1);
KurtPopulation=kurtosis(X_Big_Yields(Population)/sigma_CA,1)-3;
% event probabilities
w1=50;
w2=100;
w3=10;
w4=1;
% New fitness with higher moments
Fitness=w1*(MeanPopulation-mean_CA).^2+w2*(StDevPopulation-sigma_CA).^2+w3*SkewPopulation.^2+w4*KurtPopulation.^2;
pdeath=0.15;
pfight=0.1;
pmut=0.35;
pcross=0.35;
genimut=5;
genicross=5;
Res=zeros(20001,1);
for i=1:20000 % This time we will look into 20000 generations
    % The strongest individual is saved at index 1 in population
    [a,b]=sort(Fitness);
    aux=Population(:,b(1));
    auxF=Fitness(b(1));
    Population(:,b(1))=Population(:,1);
    Fitness(b(1))=Fitness(1);
    Population(:,1)=aux;
    Fitness(1)=auxF;
    % We will extract randselector (the individual that will be touched by
    % the event and randgenselector, the gene of the individual that will
    % be changed in the event. It is (almost) like multiplying event
    % probabilities by 1/500 (step (a) of alg. in the paper)
    randselector=rand();
    randgenselector=ceil(rand()*499)+1; % step (a)
    % The selected individual will die (step (b) of alg. in the paper) and
    % will be replaced by the strongest
    if randselector<pdeath
        Population(:,randgenselector)=Population(:,1);
    elseif randselector<pdeath+pfight
    % The selected individual will fight against another and the stronger
    % individual will replace the weaker (step (c) of alg.)
        target=ceil(rand()*500);
        winner=randgenselector;
        if Fitness(randgenselector)<Fitness(target)
            winner=target;
        end
        winnerG=Population(:,winner);
        Population(:,randgenselector)=winnerG;
        Population(:,target)=winnerG;
    % The selected individual will mutate (step (d) of alg.)
    elseif randselector<pdeath+pfight+pmut % Mutazione, step (d)
        % Mutation will touch genimut genes
        for j=1:genimut
            % A mutation is good if it does not create a repetition in the
            % individual (the realization in the pool that is selected by
            % the mutation was not present in the original individual)
            goodmut=false;
            while not(goodmut)
                target=ceil(rand()*10000); % valore del gene mutato
                goodmut=not(ismember(target,Population(:,randgenselector))); % check ammissibilit� mutazione
            end
            % Targer is a new element
            Population(ceil(rand()*100),randgenselector)=target; % la mutazione � eseguita quando ammissibile
        end
    % Last option, the selected individual will go into a crossover, step (e)
    else
        % We need a second individual for the crossover, and we make sure
        % it is different from the first one. We also make sure we are not
        % changing the strongest genoma.
        goodtarget=false;
        while not(goodtarget)
            targeti=ceil(rand()*499)+1;
            goodtarget=(targeti~=randgenselector);
        end
        % we will exchange genicross genes
        for j=1:genicross
            goodcross=false;
            % a gene is good for the exchange if it does not create repetitions.
            % This happens if: 1-the two genes point to the same
            % realization in the pool or 2-the realization in the pool that is
            % pointed to by the gene of intividual A is not already in
            % individual B genoma and vice versa
            while not(goodcross)
                targetg=ceil(rand()*100);
                if Population(targetg,targeti)==Population(targetg,randgenselector)
                    goodcross=true;
                else
                    if not(ismember(Population(targetg,targeti),Population(:,randgenselector))) && ...
                            not(ismember(Population(targetg,randgenselector),Population(:,targeti)))
                        aux=Population(targetg,targeti);
                        Population(targetg,targeti)=Population(targetg,randgenselector);
                        Population(targetg,randgenselector)=aux;
                        goodcross=true;
                    end
                end
            end
        end
    end
    % Updating fitness
    MeanPopulation=mean(X_Big_CashAccount(Population));
    StDevPopulation=std(X_Big_Yields(Population),1);
    SkewPopulation=skewness(X_Big_Yields(Population)/sigma_CA,1);
    KurtPopulation=kurtosis(X_Big_Yields(Population)/sigma_CA,1)-3;
    Fitness=w1*(MeanPopulation-mean_CA).^2+w2*(StDevPopulation-sigma_CA).^2+w3*SkewPopulation.^2+w4*KurtPopulation.^2;
    if (Fitness(1)~=Res(i))
        disp(Fitness(1));
    end
    if (mod(i,100)==0)
        disp(i)
    end
    Res(i+1)=Fitness(1);
end
% As usual, a couple of plots to look into population
subplot(1,2,1)
normplot(log(X_Big_CashAccount(Population(:,1))))
subplot(1,2,2)
histfit(log(X_Big_CashAccount(Population(:,1))),10)
figure()
plot(X_Big(Population(:,1),:)')
disp(mean(X_Big_CashAccount(Population(:,1))))
disp(std(log(X_Big_CashAccount(Population(:,1))),1))
disp(skewness(X_Big_Yields(Population(:,1))/sigma_CA,1))
disp(kurtosis(X_Big_Yields(Population(:,1))/sigma_CA,1))
