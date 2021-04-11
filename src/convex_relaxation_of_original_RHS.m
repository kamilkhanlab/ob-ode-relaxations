%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function evaluates McCormick convex relaxations for the original    %
%RHS function f_i(t, p, x) (see Assumption 3 in the article)              %
%                                                                         %
%Last modified by Yingkai Song 08/18/2020                                 %                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcvi = convex_relaxation_of_original_RHS(t,p,x,xL,xU,pL,pU,i,original_RHS)
    %retrieve the number of state variables and uncertain parameters
    nx = length(x);
    np = length(p);
    %preallocate vectors of McCormick objects
    xMC(1:nx) = McCormick(0,0,0,0);
    pMC(1:np) = McCormick(0,0,0,0);
    %construct McCormick objects xMC of (x, xL, xU) and pMC of (p, pL, pU)
    for j = 1:nx
        xMC(j) = McCormick(xL(j),xU(j),x(j),x(j));
    end
    for j = 1:np
        pMC(j) = McCormick(pL(j),pU(j),p(j),p(j));
    end
    %apply McCormick relaxation method to the original RHS and return the
    %convex relaxation
    fiMC = original_RHS(t,pMC,xMC,i);
    fcvi = fiMC.convex;
end