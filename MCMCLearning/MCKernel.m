function g = MCKernel(dt, para)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the value of kernel function at different time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dt = t_current - t_hist(:);
distance = dt;

switch para.kernel
    case 'exp'
        g = para.w * exp(-para.w * distance);
        g(g>1) = 0;
        
    case 'gauss'
        g = exp(-(distance.^2)./(2*para.w^2))./(sqrt(2*pi)*para.w);
        
    otherwise
        disp('Error: please assign a kernel function!');
end
    
