% Example 2: we will look in vectors of random variable that are used
% to simulate a geometric Brownian motion. The Browian motion can be
% seen as an asset, yield 4%, volatility 15%. The timespan of the
% simulation is 1 year (twelve monthly realization). We have a 10000 paths
% pool and will select a subsample of 100 paths in the pool to obtain
% the same performance of the original pool
rng (13,'twister')
% Browinan motion parameters
Nper=12;
sigma=0.15;
r=0.04;
dt=1/Nper;
% Creation of the Pool
X_Big_Processes=r*dt+randn(10000,Nper)*sigma*(dt^0.5)-0.5*(sigma^2)*dt;
X_Big_Yields=sum(X_Big_Processes,2);
X_Big_Asset=exp(X_Big_Yields);
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
% Initialize fitness
MeansOfPopulation=mean(X_Big_Asset(Population));
StDevsOfPopulation=std(X_Big_Yields(Population));
w1=10;    
w2=1;
Fitness=w1*(MeansOfPopulation-exp(r)).^2+w2*(StDevsOfPopulation-sigma).^2;
% Probability of events
pdeath=0.1;
pfight=0.1;
pmut=0.4;
pcross=0.3;
Res=zeros(10001,1);
genimut=5;
genicross=5;
for i=1:10000 % we will look into 10000 generations
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
    randgenselector=ceil(rand()*499)+1;
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
            % Il target è un nuovo elemento
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
    % Upating fitness
    MeansOfPopulation=mean(X_Big_Asset(Population));
    StDevsOfPopulation=std(X_Big_Yields(Population));
    Fitness=w1*(MeansOfPopulation-exp(r)).^2+w2*(StDevsOfPopulation-sigma).^2;
    if (Fitness(1)~=Res(i))
        disp(Fitness(1));
    end
    if (mod(i,100)==0)
        disp(i)
    end
    Res(i+1)=Fitness(1);
end
% A couple of plots to look into population
subplot(1,2,1)
normplot(log(X_Big_Asset(Population(:,1))))
subplot(1,2,2)
histfit(log(X_Big_Asset(Population(:,1))),15)
figure()
plot(cumsum([zeros(100,1),X_Big_Processes(Population(:,1),:)],2)')
disp(mean(X_Big_Asset(Population(:,1))))
disp(std(log(X_Big_Asset(Population(:,1)))))
