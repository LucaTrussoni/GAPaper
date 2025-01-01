% Example 2+: same as example 2 (GA2b.m) but we introduce
% 3rd and 4th moment in the fitness. We will need more epochs
rng (3,'twister')
Nper=12;
sigma=0.15;
r=0.04;
dt=1/Nper;
X_Big_Processes=r*dt+randn(10000,Nper)*sigma*(dt^0.5)-0.5*(sigma^2)*dt;
% twelve monthly returns, tannual rate 4%, volatility 15%
X_Big_Yields=sum(X_Big_Processes,2);
X_Big_CashAccount=exp(X_Big_Yields);
Population=ceil(rand(100,500)*10000);
% make sure all vectors are individuals
for i=1:500
    while length(unique(Population(:,i)))<100
        Population(:,i)=ceil(rand(100,1)*10000);
    end
end
% 1000 individuals with 100 genes each
MeanPopulation=mean(X_Big_CashAccount(Population));
StDevPopulation=std(X_Big_Yields(Population),1);
SkewPopulation=skewness(X_Big_Yields(Population)/sigma,1);
KurtPopulation=kurtosis(X_Big_Yields(Population)/sigma,1)-3;
w1=125;
w2=100;
w3=10;
w4=1;
Fitness=w1*(MeanPopulation-exp(r)).^2+w2*(StDevPopulation-sigma).^2+w3*SkewPopulation.^2+w4*KurtPopulation.^2;
pdeath=0.1;
pfight=0.1;
pmut=0.4;
pcross=0.3;
Res=zeros(10001,1);
genimut=5;
genicross=5;
for i=1:20000 % epochs
    % Save the best individual in the first place
    [a,b]=sort(Fitness);
    aux=Population(:,b(1));
    auxF=Fitness(b(1));
    Population(:,b(1))=Population(:,1);
    Fitness(b(1))=Fitness(1);
    Population(:,1)=aux;
    Fitness(1)=auxF;
    % We will extract randselector (the individual that will be touched by
    % the event and randgenselector, the gene of the individual that will
    % be changed in the event.
    randselector=rand();
    randgenselector=ceil(rand()*499)+1;
    % The selected individual will die and
    % will be replaced by the strongest
    if randselector<pdeath
        Population(:,randgenselector)=Population(:,1);
    % The selected individual will fight against another and the stronger
    % individual will replace the weaker
    elseif randselector<pdeath+pfight
        target=ceil(rand()*500);
        winner=randgenselector;
        if Fitness(randgenselector)<Fitness(target)
            winner=target;
        end
        winnerG=Population(:,winner);
        Population(:,randgenselector)=winnerG;
        Population(:,target)=winnerG;
    % The selected individual will mutate
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
            % The target is a new element
            Population(ceil(rand()*100),randgenselector)=target;
        end
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
    SkewPopulation=skewness(X_Big_Yields(Population)/sigma,1);
    KurtPopulation=kurtosis(X_Big_Yields(Population)/sigma,1)-3;
    % New fitness with highest moments
    Fitness=w1*(MeanPopulation-exp(r)).^2+w2*(StDevPopulation-sigma).^2+w3*SkewPopulation.^2+w4*KurtPopulation.^2;
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
normplot(log(X_Big_CashAccount(Population(:,1))))
subplot(1,2,2)
histfit(log(X_Big_CashAccount(Population(:,1))),10)
figure()
plot(cumsum([zeros(100,1),X_Big_Processes(Population(:,1),:)],2)')
disp(mean(X_Big_CashAccount(Population(:,1))))
disp(std(log(X_Big_CashAccount(Population(:,1))),1))
disp(skewness(log(X_Big_CashAccount(Population(:,1))),1))
disp(kurtosis(log(X_Big_CashAccount(Population(:,1))),1))
