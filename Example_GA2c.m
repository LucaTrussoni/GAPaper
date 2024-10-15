% Example 2: we will look in vectors of random variable that are used
% to simulate a an Orstein-Uhlembeck process (the random part of a Vasicek
% IR model). The process has equation dx=-a x dt + sigma dW. Again, we will
% simulate the process for 1 year (twelve monthly steps).We have a 10000
% simulated path pool and will select a subsample of 100 paths in the pool to obtain
% the same performance of the original pool
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
X_Big_Rendimenti=log(X_Big_CashAccount);
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
MeansOfPopulation=mean(X_Big_CashAccount(Population));
StDevsOfPopulation=std(X_Big_Rendimenti(Population),1);
w1=100;
w2=10;
Fitness=w1*(MeansOfPopulation-mean_CA).^2+w2*(StDevsOfPopulation-sigma_CA).^2;
% event probabilities
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
    % The selected individual will fight against another and the stronger
    % individual will replace the weaker (step (c) of alg.)
    elseif randselector<pdeath+pfight
        target=ceil(rand()*500);
        winner=randgenselector;
        if Fitness(randgenselector)<Fitness(target)
            winner=target;
        end
        winnerG=Population(:,winner);
        Population(:,randgenselector)=winnerG;
        Population(:,target)=winnerG;
    % The selected individual will mutate (step (d) of alg.)
    elseif randselector<pdeath+pfight+pmut
        % Mutation will touch genimut genes
        for j=1:genimut
            % A mutation is good if it does not create a repetition in the
            % individual (the realization in the pool that is selected by
            % the mutation was not present in the original individual)
            goodmut=false;
            while not(goodmut)
                target=ceil(rand()*10000);
                goodmut=not(ismember(target,Population(:,randgenselector)));
            end
            Population(ceil(rand()*100),randgenselector)=target;
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
    MeansOfPopulation=mean(X_Big_CashAccount(Population));
    StDevsOfPopulation=std(X_Big_Rendimenti(Population),1);
    Fitness=w1*(MeansOfPopulation-mean_CA).^2+w2*(StDevsOfPopulation-sigma_CA).^2;
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
normplot(X_Big_Rendimenti(Population(:,1)))
subplot(1,2,2)
histfit(X_Big_Rendimenti(Population(:,1)),10)
figure()
plot(X_Big(Population(:,1),:)')
disp(mean(X_Big_CashAccount(Population(:,1))))
disp(std(X_Big_Rendimenti(Population(:,1)),1))
