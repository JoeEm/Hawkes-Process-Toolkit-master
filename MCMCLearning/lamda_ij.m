%calculating the value of lamda_ij(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get rid of the m-basis function
%model.a -> a(t;parameters)
%model.mu -> Gaussian distribution with paramenters(mean,sig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = lamda_ij(Seqs,model,i_son,j_parent,timestamp,alg)

% initial 
muest = model.mu;
a = model.A;

muest = muest/sum(muest);

%D = size(Aest, 1);

tic;
        
Loglike = 0; % negative log-likelihood
lambdai = 0;


% E-step: evaluate the responsibility using the current parameters    
for c = 1:length(Seqs)
    Time = Seqs(c).Time;
    Event = Seqs(c).Mark;
    Tstart = Seqs(c).Start;

    if isempty(alg.Tmax)
        Tstop = Seqs(c).Stop;
    else
        Tstop = alg.Tmax;
        indt = Time < alg.Tmax;
        Time = Time(indt);
        Event = Event(indt);
    end

    dT = Tstop - Time;
    %modified section
    GK = Kernel_Integration(dT, model);
    GK = mean(GK);

    Nc = length(Time);

    for i = 1:Nc
        ui = Event(i);
        %verified section
        if (Event(i) ~= i_son)
            %lambdai = lambdai + 0.01;
            continue; 
        end
        ti = Time(i);             

        lambdai = muest(ui);

        if i>1
            
            %get rid of the cause event other than j.
            index = find(Event == j_parent);
            for k = 1:length(index)
                tj = Time(index(k));
                if (ti <= tj)
                   continue; 
                end

                amplifier = 10000;
                dt = ti - tj;
                gij = MCKernel(dt, model);
                %additonal part
                gij = amplifier*gij*a(i_son,j_parent,timestamp);  %sparsity problem
                %gij = gij*a(i_son,j_parent,timestamp); 

                lambdai = lambdai + sum(gij);
            end
        end
    end
        
end
 
value = lambdai;
if value <0
    value = 0;
end

end