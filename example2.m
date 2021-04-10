%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This example is modified from an anaerobic digestion model by Bernard    %
%et al. (2001)                                                            %
%                                                                         %
%Last modified by Yingkai Song 08/18/2020                                 % 
%                                                                         %
%Reference                                                                %
%         Bernard, O., Hadj-Sadok, Z., Dochain, D., Genovesi, A.,Steyer,  % 
%         J.P.: Dynamical model development and parameter identification  %
%         for an anaerobic wastewater treatment process. Biotechnology and% 
%         Bioengineering 75(4), 424-438 (2001)                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter bounds, values, tspan and other settings 

%bounds of uncertain parameters
pL = [22.14,80,238,30.6,313.6,423,4.28,0.84];
pU = [62.14,146.5,298,70.6,373.6,483,14.28,1.16];

%parameter values of interest
p = [42.14,116.5,268,50.6,343.6,450,9.28,1.0];

%time horizon of the original parametric ODE system
tspan = [0 2];

%settings for ODE solvers and NLP solver fmincon
ODE_solver = @ode23;
ODE_solver_options = odeset('RelTol',1e-4,'AbsTol',1e-4);
fmincon_options = optimoptions('fmincon','Algorithm','interior-point','OptimalityTolerance',1e-6,'Display','off');

%% Compute new optimization-based state relaxations
[t,xAug] = compute_state_relaxations(p,pL,pU,tspan,@original_initial_value,@original_RHS,ODE_solver,ODE_solver_options,fmincon_options);
%%note: this run takes ~350 seconds on a desktop with two 3.00 GHz Intel Core i7-9700 CPUs, 16.0 GB of RAM 
%% Original initial-condition and RHS functions

%initial-condition function x_0(p)
function x0 = original_initial_value(p)
    nx = 5;                  %number of state variables
    x0 = repmat(p(1),nx,1);  %always use ``repmat" to preallocate the output vector
    x0(1) = 0.5+0.*p(1);     %always add trivial p-dependence if x0 is a constant 
    x0(2) = p(8);
    x0(3) = 1+0.*p(1);
    x0(4) = 5+0.*p(1);
    x0(5) = 40+0.*p(1);
end

%RHS function of the original parametric ODE system of the form
%f_i(t,p,x), for each component i
function fi = original_RHS(~,p,x,i)
    mu1 = 1.2.*(1./(1+7.1.*(1./x(3))));
    mu2 = 0.74.*(1./(1+p(7).*(1./x(4))+x(4).*(1./256)));
    phi = x(5)+x(4)-34+p(6).*(1./19.8).*mu2.*x(2);
    q = 19.8.*(x(5)+x(4)-50-0.5.*phi);
    if i == 1
        fi = (mu1-0.2).*x(1);
    elseif i == 2
        fi = (mu2-0.2).*x(2);
    elseif i == 3
        fi = 0.4.*(5-x(3))- p(1).*mu1.*x(1);
    elseif i == 4
        fi = 0.4.*(80-x(4))+p(2).*mu1.*x(1)-p(3).*mu2.*x(2);
    elseif i == 5
        fi = -0.4.*x(5)-q+p(4).*mu1.*x(1)+p(5).*mu2.*x(2);
    end
            
end