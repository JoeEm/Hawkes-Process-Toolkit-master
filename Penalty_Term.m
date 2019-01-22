%Three Types of regularizers:
%(1)Local Indepence
%(2)Temporal Sparsity
%(3)Pairwise Smilarity

function PenaltyValue = Penalty_Term(model,alg)
    A = model.A;  %3-D martix
    %refer to the experiment of the paper
    As = 10;
    Ag = 100;
    Ap = 1000;
    
    
    R1 = 0;
    R2 = 0;
    R3 = 0;
    
    if (alg.Sparse==1)
        for i = 1:size(A,1)
            for j = 1:size(A,1)
                R1 = R1 + norm(A(i,:,j));  %2-norm regularizer
                R2 = R2 + sum(abs(A(i,:,j)));
            end 
        end
        R1 = As*R1;
        R2 = Ag*R2;
    end
    
    if (alg.Sparse==1 && alg.GroupSparse ==1) 
       for i = 1:size(A,1) 
            for j = 1:size(A,1)
                
            end 
        end
    end
    %return the penalty-term value
    PenaltyValue = As*R1 + Ag*R2 + Ap*R3;
end