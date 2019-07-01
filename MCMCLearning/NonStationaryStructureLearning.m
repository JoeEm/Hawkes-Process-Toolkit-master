%using MCMC algorithnm to estimate the non-stationary a_ij(t;parameters)
%options.N = 500; % the number of sequences                       2000
%options.Nmax = 1000; % the maximum number of events per sequence 500
%options.Tmax = 400; % the maximum size of time window            200
%In the another half of the Seqs, it generates data
%corresponding to a different structure patctern.
%clear the workplace
% N = length of each seqs.
%
clc
clear; 

startpoint = 1;

%Intialization of parameter
%for parameters 'alpha' from f() alpha -> exp(k(i))
% k(t) = c*exp( (-1)*d(t-u) )
alpha(1)=10;
final_alpha(startpoint) = alpha(1);
sig_alpha = 0.01;

%for parameters 'phai' g()  phai -> normrnd(phai(i),sig)
% k(t) = c*exp( (-1)*d(t-u) )
phai(1)=5;
final_phai(startpoint) = phai(1);
sig_phai = 0.01;
%for parameters from decay-kernel  w -> normrnd(w(i),sig)
w(1)=10;
final_w(startpoint) = w(1);
sig_w = 0.001;
%hyperparamenter weight of society activity
%weight = 2;

%function f(EventCounter,TimelengthofSeqs,alpha)
f = @(EventCounter,TimelengthofSeqs,alpha) alpha*EventCounter/TimelengthofSeqs; 

%function g(lamda_ij(t); phai) signoid function
g = @(lamda_ij, phai) 1/(1+exp(-(lamda_ij+phai))); 

%preparing simulation data
%[Seqs,alg] = GeneratingSimulationData();
load NonSatationaryTestSeqs.mat
[NewSeqs1, NewSeqs2]  = BreakSeqs(Seqs1);
%Seqs = Seqs1;
Seqs = Seqs1;
alg = options;

%Initializing the a_ij(t;parameters)
D = 8;
Tmax = options.Tmax/2;  %verified section
A = rand(D,D,startpoint);

%for parameters mu as a exogenous base intensity   mu -> normrnd(mu(i),sig)
mu=rand(D,1)+10;
final_mu(:,startpoint) = mu(:,1);
sig_mu = 0.01; 

counter = 1;
GlobalFirstTimeFlag = 0;

for t = startpoint:Tmax/10  %add the interval
    GlobalFirstTimeFlag = 1;
    for seqNum = 1:size(Seqs,2)
        index = find( (Seqs(seqNum).Time<t*Tmax/10) );
        if (length(index) < 1)
            continue;
        end
        for i=1:length(index)
            u=rand;
            if (GlobalFirstTimeFlag == 1)
                alpha(i) = final_alpha(t);
                phai(i) = final_phai(t);
                mu(:,i) = final_mu(t);
                w(i) = final_w(t);
                GlobalFirstTimeFlag=0;
            end
            New_alpha = normrnd(alpha(i),sig_alpha);
            New_phai = normrnd(phai(i),sig_phai);
            New_mu = normrnd(mu(:,i),sig_mu);
            New_w = normrnd(w(i),sig_w);
            
            Pre_model.A = A;
            Pre_model.mu = mu(:,i);                                                                                                                                                                                                                                                                                                               
            Pre_model.landmark = 0;
            Pre_model.kernel = 'gauss';
            Pre_model.w = w(i);
             
            %function f(EventCounter,TimelengthofSeqs,alpha)
            %calculating previous-likelihood
            NewSeqs = CutSeqs(Seqs(seqNum), length(index));

            Pre_Loglike = Loglike_NonSataionaryA( NewSeqs, Pre_model, alg);
               
            Func_f = f( size(Seqs(seqNum).Mark ,2), Seqs(seqNum).Stop-Seqs(seqNum).Start, alpha(i));
            tmp = 0;
            for child_node = 1:D
                for  parent_node = 1:D
                    %calculating previous-likelihood
                    %function g(lamda_ij(t); phai) signoid function
                    lamda_Value(child_node,parent_node) = lamda_ij(NewSeqs,Pre_model,child_node,parent_node,t,alg);
                    tmp = tmp + lamda_Value(child_node,parent_node);
                    %Func_g(child_node,parent_node) = g(lamda_Value(child_node,parent_node),phai(i)); 
                end          
            end
            
            for child_node = 1:D
                for  parent_node = 1:D
                    if (rand < lamda_Value(child_node,parent_node)/sum(lamda_Value(child_node,:)))
                        %New_A(child_node,parent_node,t) = lamda_Value(child_node,parent_node)/tmp;  %verified section
                        %New_A(child_node,parent_node,t) = lamda_Value(child_node,parent_node);
                        New_A(child_node,parent_node,t) = 1;
                    else
                        New_A(child_node,parent_node,t) = 0;   %
                    end
                end          
            end
            
            
            
            %calculating the Cur_Loglikehood
            Cur_model.A = New_A;
            Cur_model.mu = New_mu;
            Cur_model.landmark = 0;
            Cur_model.kernel = 'gauss';
            Cur_model.w = New_w;
            
            Cur_Loglike = Loglike_NonSataionaryA( NewSeqs, Cur_model,alg);
            
            New_alpha = normrnd(alpha(i),sig_alpha);
            New_phai = normrnd(phai(i),sig_phai);
            New_mu = normrnd(mu(:,i),sig_mu);
            New_w = normrnd(w(i),sig_w);
            
            preToCur_alpha = normpdf(New_alpha,alpha(i),sig_alpha);
            preToCur_phai = normpdf(New_phai,phai(i),sig_phai);
            preToCur_mu = mean(normpdf(New_mu,mu(:,i),sig_mu));
            preToCur_w = normpdf(New_w,w(i),sig_w);
            preToCurValue = preToCur_alpha*preToCur_phai*preToCur_mu*preToCur_w;
            
            curToPre_alpha = normpdf(alpha(i),New_alpha,sig_alpha);
            curToPre_phai = normpdf(phai(i),New_phai,sig_phai);
            curToPre_mu = mean(normpdf(mu(:,i),New_mu,sig_mu));
            curToPre_w = normpdf(w(i),New_w,sig_w);
            curToPreValue = curToPre_alpha*curToPre_phai*curToPre_mu*curToPre_w;
            
            if u<min(1,(Cur_Loglike*curToPreValue)/(Pre_Loglike*preToCurValue))  
                alpha(i+1) = New_alpha;
                phai(i+1) = New_phai;
                mu(:,i+1) = New_mu;
                w(i+1) = New_w;
            else
                alpha(i+1) = alpha(i);
                phai(i+1) = phai(i);  
                mu(:,i+1) = mu(:,i);
                w(i+1) = w(i); 
            end
        end
        if (length(alpha) ~= 1)
            tmp_Alpha = alpha(i+1);
            tmp_Phai = phai(i+1);
            tmp_Mu = mu(:,i+1);
            tmp_W = w(i+1);
            alpha= 0;
            phai = 0;
            mu = 0;
            w = 0;
            alpha= tmp_Alpha;
            phai = tmp_Phai;
            mu = tmp_Mu;
            w= tmp_W;
        else
            tmp_Alpha = alpha(i);
            tmp_Phai = phai(i);
            tmp_Mu = mu(:,i);
            tmp_W = w(i);
            alpha= 0;
            phai = 0;
            mu = 0;
            w = 0;
            alpha= tmp_Alpha;
            phai = tmp_Phai;
            mu = tmp_Mu;
            w= tmp_W;
        end
    end
    A(:,:,t+1) = New_A(:,:,size(New_A,3));
    
    if (length(alpha) ~= 1)
        final_alpha(t+1) = alpha(i+1);
        final_phai(t+1) = phai(i+1);
        final_mu(:,t+1) = mu(:,i+1);
        final_w(t+1) = w(i+1);
    else
        final_alpha(t+1) = alpha(1);
        final_phai(t+1) = phai(1);
        final_mu(:,t+1) = mu(:,1);
        final_w(t+1) = w(1);
    end
    dis_alpha = final_alpha(length(final_alpha))
    dis_phai = final_phai(length(final_phai))
    dis_mu =final_mu(:,size(final_mu,2))
    dis_w = final_w(length(final_w))
    dis_mean_A = mean(sum(sum(sum(New_A(:,:,t)))))
end


