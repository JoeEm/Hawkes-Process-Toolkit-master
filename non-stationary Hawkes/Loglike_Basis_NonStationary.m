function Loglike = Loglike_Basis_NonStationary( Seqs, model, alg, start, stop)
                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Hawkes processes via maximum likelihood estimation
% Different regularizers (low-rank, sparse, group sparse) of parameters and
% their combinations are considered, which are solved via ADMM.
%
% Note: do some modification on the original version to adatpt to the
% non-stationary situation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial 
Aest = model.A; 
muest = model.mu;



%D = size(Aest, 1);

tic;

        
Loglike = 0; % negative log-likelihood



% E-step: evaluate the responsibility using the current parameters    
for c = 1:length(Seqs)
    %added part comparing to the original version
    index = find(Seqs(c).Mark == 0);
    for countnum = 1:length(index)
        Seqs(c).Mark(index(countnum)) = 1;
    end
    
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
    
    
    %modified section
    %Tstart = start;
    %Tstop = stop;
    
    
    %Amu = Amu + Tstop - Tstart;
    %dbstop in Loglike_Basis_NonStationary at 52
    dT = Tstop - Time;
    GK = Kernel_Integration(dT, model);

    Nc = length(Time);

    for i = 1:Nc

        ui = Event(i);


        ti = Time(i);             

        lambdai = muest(ui);
        %pii = muest(ui);
        %pij = [];


        if i>1

            tj = Time(1:i-1);
            uj = Event(1:i-1);

            dt = ti - tj;
            gij = Kernel(dt, model);
            auiuj = Aest(uj, :, ui);
            pij = auiuj .* gij;
            lambdai = lambdai + sum(pij(:));
        end

        Loglike = Loglike - log(lambdai);

    end
    %dbstop in Loglike_Basis_NonStationary at 85
    Loglike = Loglike + (Tstop-Tstart).*sum(muest);
    Loglike = Loglike + sum( sum( GK.*sum(Aest(Event,:,:),3) ) );
    


end

%PenaltyValue = Penalty_Term(model,alg);  %add penalty term
%Loglike = Loglike + PenaltyValue;

Loglike = -Loglike;
                
        
