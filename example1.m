%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This example is a bioreaction model from Bastin and Dochain (1990)       %
%                                                                         %
%Last modified by Yingkai Song, 08/18/2020                                %
%                                                                         %
%Reference                                                                %
%         Bastin, G., Dochain, D.: On-line Estimation and Adaptive Control% 
%         of Bioreactors. Elsevier, Amsterdam (1990)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter bounds, values, tspan and other settings 

%bounds of uncertain parameters
pL = [0.4,1.2];
pU = [0.6,1.5];

%a certain p where state relaxations will be computed
p = [0.5,1.3];

%time horizon of the parametric ODE system
tspan = [0,5];

%settings for ODE solvers and NLP solver fmincon
ODE_solver = @ode23;
ODE_solver_options = odeset('RelTol',1e-4,'AbsTol',1e-4);
fmincon_options = optimoptions('fmincon','Algorithm','interior-point','OptimalityTolerance',1e-6,'Display','off');

%% Compute new optimization-based state relaxations
[t,xAug] = compute_state_relaxations(p,pL,pU,tspan,@original_initial_value,@original_RHS,ODE_solver,ODE_solver_options,fmincon_options);
%note: this run takes ~8 seconds on a desktop with two 3.00 GHz Intel Core i7-9700 CPUs and 16.0 GB of RAM 
%% Original initial-condition and RHS functions

%initial-condition function x_0(p)
function x0 = original_initial_value(p)
    nx = 2;                     %number of state variables
    x0 = repmat(p(1),nx,1);     %always use ``repmat" to preallocate the output vector
    x0(1) = 0.82+0.*p(1);       %always add trivial p-dependence if x0 is a constant 
    x0(2) = 0.8+0.*p(1);
end

%RHS function of the original parametric ODE system of the form
%f_i(t,p,x), for each component i
function fi = original_RHS(~,p,x,i)
    alpha = 0.5;
    k = 10.53;
    D = 0.36;
    Si = 5.7;
    Ks = 7.10;
    mu = p(2).*(1./(Ks.*(1./x(2))+1+p(1).*x(2)));
    if i == 1
        fi = (mu-alpha.*D).*x(1);
    else
        if i == 2
            fi = D.*(Si-x(2))-k.*mu.*x(1);
        end
    end
            
end