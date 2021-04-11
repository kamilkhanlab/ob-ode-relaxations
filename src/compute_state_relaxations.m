%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main function for computing new optimization-based relaxations for a     % 
%nonconvex parametric ODE system                                          %
%                                                                         %
%Last modified by Yingkai Song 08/18/2020                                 %
%                                                                         %
%inputs:                                                                  %
%                [pL,pU] - bounds of the uncertain parameters             %
%                      p - parameter values of interest                   %
%                  tspan - time horizon                                   %
% original_initial_value - initial-value function of the original         %
%                          parametrc ODE sytem                            %
%           original_RHS - RHS function of the original parametric ODE    %
%                          system                                         %
%             ODE_solver - solvers employed for solving the new auxiliary % 
%                          ODE system                                     %
%     ODE_solver_options - MATLAB options for the ODE solver              %
%        fmincon_options - MATLAB options for the NLP solver fmincon      %
%                                                                         %
%outputs:                                                                 %
%                      t - a vector of time steps employed by ODE solvers %
%                   xAug - a matrix that stores state variables           %
%                          [xL, xU, xcv, xcc] at each time step           % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,xAug] = compute_state_relaxations(p,pL,pU,tspan,original_initial_value,original_RHS,ODE_solver,ODE_solver_options,fmincon_options)
         %construct Interval and McCormick objects of the uncertain parameters 
         np = length(p);
         pI(1:np) = Interval(0,0);
         pMC(1:np) = McCormick(0,0,0,0);
         for i = 1:np
             pI(i) = Interval(pL(i),pU(i));
             pMC(i) = McCormick(pL(i),pU(i),p(i),p(i));
         end
         %compute the initial values of state variables [xL, xU, xcv, xcc]
         %of the auxiliary state relaxation system
         x0I = original_initial_value(pI);
         x0MC = original_initial_value(pMC);
         nx = length(x0I);
         x0L = zeros(nx,1);
         x0U = zeros(nx,1);
         x0cv = zeros(nx,1);
         x0cc = zeros(nx,1);
         for i = 1:nx
             x0L(i) = x0I(i).lower;
             x0U(i) = x0I(i).upper;
             x0cv(i) = x0MC(i).convex;
             x0cc(i) = x0MC(i).concave;
         end
         x0Aug = [x0L,x0U,x0cv,x0cc];
         %solve the new auxiliary ODE system for new optimization-based 
         %state relaxations
         [t,xAug] = ODE_solver(@(t,xAug) optimization_based_ODE_RHS(t,xAug,p,pL,pU,original_RHS,fmincon_options), tspan, x0Aug, ODE_solver_options);         
end