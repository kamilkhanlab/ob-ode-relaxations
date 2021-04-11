%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that constructs the RHS of the new state relaxation system,     %
%where the Harrison's state bounds (1977) are integrated simultanouesly,  %
%and the time derivatives of state relaxations are computed via solving   %
%convex optimization problems                                             %
%                                                                         %
%Last modified by Yingkai Song 09/20/2020                                 %
%                                                                         % 
%inputs:                                                                  %
%               t - current time step                                     %
%            xAug - a vector of state variables [xL, xU, xcv, xcc] at     %
%                   current t                                             %
%        [pL, pU] - predefined bounds of uncertain parameters             %
%               p - parameter values of interest                          %
%    original_RHS - the RHS functions of the original parametric ODE      % 
%                   system                                                %
% fmincon_options - MATLAB options for the NLP solver fmincon             %
%                                                                         %
%outputs:                                                                 %
%         dxAugdt - the time derivatives of [xL, xU, xcv, xcc]            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxAugdt = optimization_based_ODE_RHS(t,xAug,p,pL,pU,original_RHS,fmincon_options)
    %retrieve number of state variables and uncertain parameters
    nx = length(xAug)/4;
    np = length(p);
    %retrieve state bounds/relaxations from the state variable vector xAug
    xL = xAug(1:nx);
    xU = xAug(nx+1:2*nx);
    xcv = xAug(2*nx+1:3*nx);
    xcc = xAug(3*nx+1:4*nx);
    %construct interval objects xI of (xL, xU) and pI of (pL, pU)
    xI(1:nx) = Interval(0,0);
    for i = 1:nx
        xI(i) = Interval(xL(i),xU(i));
    end
    pI(1:np) = Interval(0,0);
    
    for i = 1:np
        pI(i) = Interval(pL(i),pU(i));
    end
    %preallocate the vectors of time derivatives of xL, xU, xcv, and xcc
    dxL = zeros(nx,1);
    dxU = zeros(nx,1);
    dxcv = zeros(nx,1);
    dxcc = zeros(nx,1);
    %compute time derivatives of state bounds/relaxations for each
    %component i
    for i = 1:nx
        %flatten the ith interval object (xiL, xiU) to (xiL, xiL)
        xI(i) = Interval(xL(i),xL(i));
        %apply natural interval extension to original RHS fi, to compute 
        %dxiL
        dxI = original_RHS(t,pI,xI,i);
        dxL(i) = dxI.lower;
        %flatten the ith interval object (xiL, xiU) to (xiU, xiU)
        xI(i) = Interval(xU(i),xU(i));
        %apply natural interval extension to original RHS fi, to compute 
        %dxiU
        dxI = original_RHS(t,pI,xI,i);
        dxU(i) = dxI.upper;
        %unflatten the ith interval object xiI
        xI(i) = Interval(xL(i),xU(i));
        
        %compute time derivatives of state relaxations
        
        %lower and upper bounds for decision variables alpha
        lb = -ones(nx,1);
        ub = ones(nx,1);
        %equality constraint: alpha_i = -1
        Aeq = zeros(1,nx);
        Aeq(i) = 1;
        beq = -1;
        %construct and solve a convex optimization problem for dxicv
        [~,dxcv(i)] = fmincon(@(alpha) convex_relaxation_of_original_RHS(t,p,linear_transformation(alpha,xcv,xcc),xL,xU,pL,pU,i,original_RHS), zeros(nx,1),[],[],Aeq,beq,lb,ub,[],fmincon_options);
        %equality constraint: alpha_i = 1
        beq = 1;
        %construct and solve a convex optimization problem for dxicc
        [~,ndxcc] = fmincon(@(alpha) -concave_relaxation_of_original_RHS(t,p,linear_transformation(alpha,xcv,xcc),xL,xU,pL,pU,i,original_RHS), zeros(nx,1),[],[],Aeq,beq,lb,ub,[],fmincon_options);
        dxcc(i) = -ndxcc;
        %discrete jump in Scott--Barton state relaxation framework (2013)
        if xcv(i)<=xL(i)
            dxcv(i) = max([dxcv(i),dxL(i)]);
        end
        if xcc(i)>=xU(i)
            dxcc(i) = min([dxcc(i),dxU(i)]);
        end  
    end
    dxAugdt = [dxL;dxU;dxcv;dxcc];
end