function [alpha1,alpha2] = fAlpha(epsilon,epsilonTarget,alpha1min,alpha2min)
global countAlpha;
if isempty(countAlpha)
    countAlpha=0;
end
alpha1 = 1;
alpha2 = 1;
if epsilon == epsilonTarget
    countAlpha = countAlpha + 1;
    if alpha1 > alpha1min
        alpha1 =  1/(1+countAlpha/20);
    else 
        alpha1 = alpha1min;
    end
    if alpha2 > alpha2min
        alpha2 =  1/(1+countAlpha/20);
    else 
        alpha2 = alpha2min;
    end
end