function Loglike = Loglike_NonSataionaryA( Seqs, model, alg )
                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get rid of the m-basis function
%model.a -> a(t;parameters)
%model.mu -> Gaussian distribution with paramenters(mean,sig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial 
muest = model.mu;
a = model.A;

tic;        
Loglike = 0; % negative log-likelihood



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
        ti = Time(i);             

        lambdai = muest(ui);
        %pii = muest(ui);
        %pij = [];

        if i>1

            tj = Time(1:i-1);
            %uj = Event(1:i-1); 

            dt = ti - tj;
            gij = Kernel(dt, model);
            %additonal part 
            ti = floor(ti/50);
            if (ti == 0)
                ti = ti + 1;
            end
            AvergeofA = sum(a(ui,:,ti))/size(a,2);  %approximate it by average value.
            gij = gij*AvergeofA;

            lambdai = lambdai + sum(gij);
        end
        
        if (lambdai>0)
            Loglike = Loglike + log(lambdai);
        end

    end

   Loglike = Loglike - (Tstop-Tstart).*sum(muest);
   %Loglike = Loglike + sum( sum( GK.*sum(Aest(Event,:,:),3) ) );
   AvergeOfA = sum(sum(sum(a)))/(size(a,1)*size(a,2)*size(a,3)); %take the average as an approximation. 

   Loglike = Loglike - GK*AvergeOfA;

end

Loglike = Loglike;

end
