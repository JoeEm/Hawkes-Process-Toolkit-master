%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function :NonStationaryDataGenerating(void)
% 1.Generating nonStationary Data...(Two patterns).
% 2.Displaying the two different patterns respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Seqs=DataGeneratingNonStationary(alg1)


options.N = 700; % the number of sequences 2000
options.Nmax = 500; % the maximum number of events per sequence
options.Tmax = 200; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 100; % the number of generations for branch processing
D = 8; % the dimension of Hawkes processes


disp('Approximate simulation of Hawkes processes via branching process')
disp('Complicated gaussian kernel')
para1.kernel = 'gauss'; % the type of kernels per impact function
para1.w = 2; % the bandwidth of gaussian kernel
para1.landmark = 0:4:12; % the central locations of kernels
L = length(para1.landmark);

% initialize ground truth parameters
para1.mu = rand(D,1)/D;
para1.A = zeros(D, D, L);
for l = 1:L
    para1.A(:,:,l) = (0.5^l)*(0.5+ones(D));
end
mask = rand(D).*double(rand(D)>0.7);
para1.A = para1.A.*repmat(mask, [1,1,L]);
para1.A = 0.25*para1.A./max(abs(eig(sum(para1.A,3)))); % ensure the stationarity of Hawkes process
tmp = para1.A;
para1.A = reshape(para1.A, [D, L, D]);
for di = 1:D
    for dj = 1:D
        phi = tmp(di, dj, :);
        para1.A(dj, :, di) = phi(:);
    end
end


% initialize ground truth parameters
para1.mu = rand(D,1)/D;
para1.A = zeros(D, D, L);
for l = 1:L
    para1.A(:,:,l) = (0.5^l)*(0.5+ones(D));
end
mask = rand(D).*double(rand(D)>0.7);
para1.A = para1.A.*repmat(mask, [1,1,L]);
para1.A = 0.25*para1.A./max(abs(eig(sum(para1.A,3)))); % ensure the stationarity of Hawkes process
tmp = para1.A;
para1.A = reshape(para1.A, [D, L, D]);
for di = 1:D
    for dj = 1:D
        phi = tmp(di, dj, :);
        para1.A(dj, :, di) = phi(:);
    end
end

para1.A = (para1.A)*10;


% two simulation methods
Seqs1 = Simulation_Branch_HP(para1, options);
model1 = Initialization_Basis(Seqs1); 
% learning the model of Seqs1 by MLE-SGL

figure
A1_T=ImpactFunc(para1, options);
subplot(221);
imagesc(A1_T);
colorbar
title('Ground truth of infectivity_A1')

model_MLE_SGL = Learning_MLE_S_Basis( Seqs1, model1, alg1);
A1 = ImpactFunc(model_MLE_SGL, options);


subplot(222);
imagesc(A1);
colorbar
title('Estimation of infectivity_A1')
%Seqs1 = Simulation_Thinning_HP(para1, options);


%generating pattern two
for di = 1:D
    for dj = 1:D
        if (...
            (di == 1 && (dj == 1 || dj ==2 || dj ==3 || dj ==4)) ||... 
            (di == 2 && (dj == 1 || dj ==2 || dj ==3 || dj ==4)) ||...
            (di == 3 && (dj == 1 || dj ==2 || dj ==3 || dj ==4)) ||...
            (di == 4 && (dj == 1 || dj ==2 || dj ==3 || dj ==4)) ||...
            (di == 5 && (dj == 5 || dj ==6 || dj ==7 || dj ==8)) ||...
            (di == 6 && (dj == 5 || dj ==6 || dj ==7 || dj ==8)) ||...
            (di == 7 && (dj == 5 || dj ==6 || dj ==7 || dj ==8)) ||...
            (di == 8 && (dj == 5 || dj ==6 || dj ==7 || dj ==8)))
                        
                para1.A(dj, :, di) = 0;
        end
    end
end

%para1.A = (para1.A)*10;

para1.w = 2; % the bandwidth of gaussian kernel
para1.landmark = 0:4:12; % the central locations of kernels


options.N = 300; % the number of sequences
options.Nmax = 500; % the maximum number of events per sequence
options.Tmax =200; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 100; % the number of generations for branch processing
D = 8; % the dimension of Hawkes processes


Seqs2 = Simulation_Branch_HP(para1, options);
model2 = Initialization_Basis(Seqs1); 

A2_T=ImpactFunc(para1, options);
subplot(223);
imagesc(A2_T);
colorbar
title('Ground truth of infectivity_A2')

% learning the model of Seqs1 by MLE-SGL
alg1.Sparse=1;
alg1.GroupSparse=1; 
model_MLE_SGL = Learning_MLE_S_Basis( Seqs2, model2, alg1 );
A2 = ImpactFunc(model_MLE_SGL, options);

subplot(224);
imagesc(A2);
colorbar
title('Estimation of infectivity_A2')
% Visualize the infectivity matrix (the adjacent matrix of Granger causality graph)
  
Seqs = [Seqs1,Seqs2]; %combining the two different patterns into one  



end
