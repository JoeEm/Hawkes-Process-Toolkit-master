% Metropolis(-Hastings) algorithm
% true (target) is the loglikelihood of lamda(t;parameters). 
% proposal (sample) pdf is Gaussian or exponential where we can sample.
%% 
clc
clear; 



X(1)=0; 
N=1e4;
p = @(x) 0.3*exp(-0.2*x.^2) + 0.7*exp(-0.2*(x-10).^2);  %Note: p must be a pdf!!!
dx=0.5; xx=-10:dx:20; fp=p(xx); plot(xx,fp) % plot the true p(x)
%% MH algorithm
sig=(10); 
for i=1:N-1
    u=rand;
    x=X(i); 
    xs=normrnd(x,sig); % new sample xs based on existing x from proposal pdf.
    pxs=p(xs);
    px=p(x); 
    qxs=normpdf(xs,x,sig); % qxs ->q(old->nex)  
    qx=normpdf(x,xs,sig); %  qx ->q(new->old)
     if u<min(1,pxs*qx/(px*qxs))  % case 1: pesudo code
%     if u<min(1,pxs/(px))        % case 2: Metropolis algorithm
%     if u<min(1,pxs/qxs/(px/qx)) % case 3: independent sampler
        X(i+1)=xs;
    else
        X(it1)=x; 
    end
end
% compare pdf of the simulation result with true pdf.
N0=1;  close all;figure; %N/5; 
nb=histc(X(N0+1:N),xx); 
bar(xx+dx/2,nb/(N-N0)/dx); % plot samples.
A=sum(fp)*dx; 
hold on; plot(xx,fp/A,'r') % compare.
% figure(2); plot(N0+1:N,X(N0+1:N)) % plot the traces of x.

% compare cdf with true cdf.
F1(1)=0;
F2(1)=0;
for i=2:length(xx) 
  F1(i)=F1(i-1)+nb(i)/(N-N0); 
  F2(i)=F2(i-1)+fp(i)*dx/A;
end

figure
plot(xx,[F1' F2'])
max(F1-F2) % this is the true possible measure of accuracy.