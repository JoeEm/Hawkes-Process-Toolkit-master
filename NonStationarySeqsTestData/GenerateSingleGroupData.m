%In the another half of the Seqs, it generates data
%corresponding to a different structure pattern.

function [Seqs1,para] = GenerateSingleGroupData(options,D,Group,amplifier,MaskValue, A, flag)


%for Tmax 1 = 1honr, 1 week = 24*7 = 168, 

disp('Approximate simulation of Hawkes processes via branching process')
disp('Complicated gaussian kernel')
para.kernel = 'gauss'; % the type of kernels per impact function
para.w = 0.2; % the bandwidth of gaussian kernel
para.landmark = 0:4:12; % the central locations of kernels
L = length(para.landmark);
para.mu = rand(D,1)/D;

if (flag ~= 1)
% initialize ground truth parameters
    para.A = zeros(D, D, L);
    for l = 1:L
        para.A(:,:,l) = (0.5^l)*(0.5+ones(D));
    end


    mask = rand(D).*double(rand(D)>MaskValue);
    % mask = rand(D).*double(rand(D)>0.7);
    para.A = para.A.*repmat(mask, [1,1,L]);
    para.A = 0.25*para.A./max(abs(eig(sum(para.A,3)))); % ensure the stationarity of Hawkes process
    tmp = para.A;
    para.A = reshape(para.A, [D, L, D]);
    for di = 1:D
        for dj = 1:D
            phi = tmp(di, dj, :);
            para.A(dj, :, di) = phi(:);
        end
    end
else
    para.A = A;
end
para.A = para.A * amplifier;


%Generating Simulation Data
Seqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Group', [], ...
              'Feature', []);
tic




for n = 1:options.N
    % the 0-th generation, simulate exogeneous events via Poisson processes
    History = Simulation_Thinning_Poisson(para.mu, 0, options.Tmax);
    %do some alteration on the Data structure(from 2d to 3d)
    New_History = zeros(3,length(size(History,2)));
    for q = 1:size(History,2) %two pattern situation...
        New_History(1,q) = History(1,q);
        New_History(2,q) = History(2,q);
        
        New_History(3,q) = Group;            
    end
    current_set = New_History;
    
    for k = 1:options.GenerationNum
        future_set = [];
        for i = 1:size(current_set, 2)
            ti = current_set(1,i);

            ui = current_set(2,i);
            t = 0;
            
            
            phi_t = ImpactFunction(ui, t, para);
            mt = sum(phi_t);
            
            while t<options.Tmax-ti                              
                s = random('exp', 1/mt);
                U = rand;

                phi_ts = ImpactFunction(ui, t+s, para);  %relate the Data Generation with A
                mts = sum(phi_ts);

                %fprintf('s=%f, v=%f\n', s, mts/mt);        
                if t+s>options.Tmax-ti || U>mts/mt
                    t = t+s;
                else
                    u = rand*mts;
                    sumIs = 0;
                    for d=1:length(phi_ts)
                        sumIs = sumIs + phi_ts(d);
                        if sumIs >= u
                            break;
                        end
                    end
                    index = d;

                    t = t+s;
                    %future_set=[future_set,[t+ti;index(1)]];
                    future_set=[future_set,[t+ti;index(1);Group]];
                end
        
                phi_t = ImpactFunction(ui, t, para);
                mt = sum(phi_t);
            end            
        end
        
        if isempty(future_set) || size(New_History, 2)>options.Nmax
            break
        else
            current_set = future_set;
            New_History = [New_History, current_set];
        end
    end
    
    [~, index] = sort(New_History(1,:), 'ascend');
    Seqs(n).Time = New_History(1,index);
    Seqs(n).Mark = New_History(2,index);
    
    Seqs(n).Group = New_History(3,index);
    
    Seqs(n).Start = 0;
    Seqs(n).Stop = options.Tmax;
    index = find(Seqs(n).Time<=options.Tmax);
    Seqs(n).Time = Seqs(n).Time(index);
    Seqs(n).Mark = Seqs(n).Mark(index);
    
    if mod(n, 10)==0 || n==options.N
        fprintf('#seq=%d/%d, #event=%d, time=%.2fsec\n', ...
            n, options.N, length(Seqs(n).Mark), toc);
    end
end

Seqs1 = Seqs;

end