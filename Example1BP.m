% Example 1: a genetic algorithm will select a subsample of 100
% realizations of a normal random variable from a pool of 10000
% realizations to obtain mean and variance similar to population
rng (6,'twister')
% Create the pool
X_Big=randn(10000,1);
% Population of 500 subsamples, each made of 100 realizations
Population=ceil(rand(100,500)*10000);
% Make sure all subsamples are well defined individuals, so no
% individual has to genes pointing to the same realization in the pool
% (no repetition allowed)
for i=1:500
    while length(unique(Population(:,i)))<100
        Population(:,i)=ceil(rand(100,1)*10000);
    end
end
% Initialize fitness
MeansOfIndividuals=mean(X_Big(Population));
StDevsOfIndividuals=std(X_Big(Population));
w1=10;
w2=10;
Fitness=w1*MeansOfIndividuals.^2+w2*(StDevsOfIndividuals-1).^2;
% Initialize event probabilities
pdeath=0.1;
pfight=0.05;
pmut=0.5;
pcross=0.35;
genimut=5;
genicross=5;
Res=zeros(20001,1);
for i=1:10000 % We will look into 10000 generations
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
            Population(ceil(rand()*100),randgenselector)=target;
        end
    % Last option, the selected individual will go into a crossover, step (e)
    else
        % We need a second individual for the crossover, and we make sure
        % it is different from the first one. We also make sure we are not
        % changin the strongest genoma.
        goodtarget=false;
        while not(goodtarget)
            targeti=ceil(rand()*499)+1;
            goodtarget=(targeti~=randgenselector);
        end
        % we will exchange genicross genes
        for j=1:genicross
            % a gene is good for the exchange if it does not create repetitions.
            % This happens if: 1-the two genes point to the same
            % realization in the pool or 2-the realization in the pool that is
            % pointed to by the gene of intividual A is not already in
            % individual B genoma and vice versa
            goodcross=false;
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
    % Update fitness
    MeansOfIndividuals=mean(X_Big(Population));
    StDevsOfIndividuals=std(X_Big(Population));
    w1=10;
    w2=10;
    Fitness=w1*MeansOfIndividuals.^2+w2*(StDevsOfIndividuals-1).^2;
    if (Fitness(1)~=Res(i))
        disp(Fitness(1));
    end
    if (mod(i,100)==0)
        disp(i)
    end
    Res(i+1)=Fitness(1);
end
% Make a couple of plot to look into population
subplot(1,2,1)
normplot(X_Big(Population(:,1)))
subplot(1,2,2)
histfit(X_Big(Population(:,1)),15)
