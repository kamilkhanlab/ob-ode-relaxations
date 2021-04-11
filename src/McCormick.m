%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generalized McCormick relaxations via operator overloading               %
%                                                                         %
%Overloaded operators: exp, log, -, +, .*, ./, x.^n                       %
%Accompanied interval bounds of factors are propagated via natural        %
%interval extension.                                                      %
%                                                                         %
%Last modified by Yingkai Song, 08/21/2020                                %
%                                                                         %
%Reference                                                                %
%       Scott, J.K., Stuber, M.D., Barton, P.I.: Generalized McCormick    %
%       relaxations. J. Glob. Optim. 51(4), 569-606 (2011)                %
%       Scott, J.K.: Reachability analysis and deterministic global       %
%       optimization of differential-algebraic systems. PhD thesis.       %
%       Massachusetts Institute of Technology (2012)                      %
%       Tsoukalas, A., Mitsos, A.: Multivariate McCormick relaxations.    %
%       J. Glob. Optim. 59(2-3), 633-662 (2014)                           %
%       Liberti, L., Pantelides, C.C.: Convex envelopes of monomials of   %
%       odd degree. J. Glob. Optim. 25(2), 157-168 (2003)                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef McCormick 
    properties  
        lower;        %lower bound
        upper;        %upper bound
        convex;       %convex relaxation
        concave;      %concave relaxation
    end
    methods    
        function z = McCormick(lower,upper,convex,concave)
            z.lower = lower;
            z.upper = upper;
            z.convex = convex;
            z.concave = concave;
        end
        %overload log
        function z = log(x)
            if ~isa(x,'McCormick') %when x is a number
                z = McCormick(log(x), log(x), log(x), log(x));
            else
                %perform McCormick relaxation rule for composite functions
                %(Scott et al. (2011))
                zL = log(x.lower);
                zU = log(x.upper);
                Fcc_max = x.upper;
                mid_cc = median([x.convex,x.concave,Fcc_max]);
                z_concave = log(mid_cc);
                Fcv_min = x.lower;
                mid_cv = median([x.convex,x.concave,Fcv_min]);
                if x.lower == x.upper
                    z_convex = z_concave;
                else
                    z_convex = (zU-zL)./(x.upper-x.lower).*mid_cv + (x.upper.*zL-x.lower.*zU)./(x.upper-x.lower);
                end
            z = McCormick(zL,zU,z_convex,z_concave);   
            end
        end
        %overload exp
        function z = exp(x)
            if (~isa(x,'McCormick'))
                 z = McCormick(exp(x), exp(x), exp(x), exp(x));
            else 
                zL = exp(x.lower);
                zU = exp(x.upper);
                Fcv_min = x.lower;
                mid_cv = median([x.convex,x.concave,Fcv_min]);
                z_convex = exp(mid_cv);
                Fcc_max = x.upper;
                mid_cc = median([x.convex,x.concave,Fcc_max]);
                if x.lower == x.upper
                    z_concave = z_convex;
                else
                    z_concave = (zU-zL)./(x.upper-x.lower).*mid_cc + (x.upper.*zL-x.lower.*zU)/(x.upper-x.lower);
                end
            z = McCormick(zL,zU,z_convex,z_concave);
            end                 
        end
         %overload .*
        function z = times(x,y)
            %turn non-McCormick object to McCormick object
            if (~isa(x,'McCormick'))
                x = McCormick(x,x,x,x);
            end
            if (~isa(y,'McCormick'))
                y = McCormick(y,y,y,y);
            end 
            zL = min([x.lower.*y.lower, x.lower.*y.upper, x.upper.*y.lower, x.upper.*y.upper]);
            zU = max([x.lower.*y.lower, x.lower.*y.upper, x.upper.*y.lower, x.upper.*y.upper]);
            
            %Scott-Barton product rule in Scott's PhD thesis (2012)
            alpha1 = min([y.lower.*x.convex, y.lower.*x.concave]);
            beta1 = min([y.upper.*x.convex, y.upper.*x.concave]);
            gamma1 = max([y.lower.*x.convex, y.lower.*x.concave]);
            delta1 = max([y.upper.*x.convex, y.upper.*x.concave]);
            alpha2 = min([x.lower.*y.convex,x.lower.*y.concave]);
            beta2 = min([x.upper.*y.convex, x.upper.*y.concave]);
            gamma2 = max([x.upper.*y.convex, x.upper.*y.concave]);
            delta2 = max([x.lower.*y.convex, x.lower.*y.concave]);
            z_convex = max([alpha1+alpha2-x.lower.*y.lower, beta1+beta2-x.upper.*y.upper, zL]);
            z_concave = min([gamma1+gamma2-x.upper.*y.lower, delta1+delta2-x.lower.*y.upper, zU]);
            
            %Tsoukalas-Mitsos product rule (2014) 
            %Users can comment the classic product rule above and uncomment
            %T-M rule below to obtain guaranteed tighter relaxations for 
            %bilinear terms.
            %{
            k = (y.lower-y.upper)/(x.upper-x.lower);
            zeta = (x.upper*y.upper-x.lower*y.lower)/(x.upper-x.lower);
            a = max([y.upper*x.convex+x.upper*median([y.convex,y.concave,k*x.convex+zeta])-x.upper*y.upper, y.lower*x.convex+x.lower*median([y.convex,y.concave,k*x.convex+zeta])-x.lower*y.lower]);
            b = max([y.upper*x.concave+x.upper*median([y.convex,y.concave,k*x.concave+zeta])-x.upper*y.upper, y.lower*x.concave+x.lower*median([y.convex,y.concave,k*x.concave+zeta])-x.lower*y.lower]);
            c = max([y.upper*median([x.convex,x.concave,(y.convex-zeta)/k])+x.upper*y.convex-x.upper*y.upper, y.lower*median([x.convex,x.concave,(y.convex-zeta)/k])+x.lower*y.convex-x.lower*y.lower]);
            d = max([y.upper*median([x.convex,x.concave,(y.concave-zeta)/k])+x.upper*y.concave-x.upper*y.upper, y.lower*median([x.convex,x.concave,(y.concave-zeta)/k])+x.lower*y.concave-x.lower*y.lower]);
            e = max([y.upper*x.convex+x.upper*y.convex-x.upper*y.upper, y.lower*x.convex+x.lower*y.convex-x.lower*y.lower]);
            f = max([y.upper*x.concave+x.upper*y.concave-x.upper*y.upper, y.lower*x.concave+x.lower*y.concave-x.lower*y.lower]);
            z_convex = min([a, b, c, d, e, f]);

            k = (y.upper-y.lower)/(x.upper-x.lower);
            zeta = (x.upper*y.lower-x.lower*y.upper)/(x.upper-x.lower);
            a = min([y.lower*x.convex+x.upper*median([y.convex,y.concave,k*x.convex+zeta])-x.upper*y.lower, y.upper*x.convex+x.lower*median([y.convex,y.concave,k*x.convex+zeta])-x.lower*y.upper]);
            b = min([y.lower*x.concave+x.upper*median([y.convex,y.concave,k*x.concave+zeta])-x.upper*y.lower, y.upper*x.concave+x.lower*median([y.convex,y.concave,k*x.concave+zeta])-x.lower*y.upper]);
            c = min([y.lower*median([x.convex,x.concave,(y.convex-zeta)/k])+x.upper*y.convex-x.upper*y.lower, y.upper*median([x.convex,x.concave,(y.convex-zeta)/k])+x.lower*y.convex-x.lower*y.upper]);
            d = min([y.lower*median([x.convex,x.concave,(y.concave-zeta)/k])+x.upper*y.concave-x.upper*y.lower, y.upper*median([x.convex,x.concave,(y.concave-zeta)/k])+x.lower*y.concave-x.lower*y.upper]);
            e = min([y.lower*x.convex+x.upper*y.concave-x.upper*y.lower, y.upper*x.convex+x.lower*y.concave-x.lower*y.upper]);
            f = min([y.lower*x.concave+x.upper*y.convex-x.upper*y.lower, y.upper*x.concave+x.lower*y.convex-x.lower*y.upper]);
            z_concave = max([a, b, c, d, e, f]);
            %}
            z = McCormick(zL,zU,z_convex,z_concave);                                                  
        end
        
        %overload x./y
        function z = rdivide(x,y)
            %this is overloaded by x.*(1./y)
             if (~isa(y,'McCormick'))
                 yrec = McCormick(1./y, 1./y, 1./y, 1./y);
            else
                yrecL = 1./(y.upper);
                yrecU = 1./(y.lower);
                yrec_cv_min = y.upper;
                yrec_cc_max = y.lower;
                mid_cv = median([y.convex,y.concave,yrec_cv_min]);
                mid_cc = median([y.convex,y.concave,yrec_cc_max]);
                if y.lower > 0
                    yrec_convex = 1./(mid_cv);
                    if y.lower == y.upper
                        yrec_concave = yrec_convex;
                    else
                        yrec_concave = (yrecL-yrecU)./(y.upper-y.lower).*mid_cc + (y.upper.*yrecU-y.lower.*yrecL)./(y.upper-y.lower);
                    end
                end
                if y.upper < 0
                    yrec_concave = 1./(mid_cc);
                    if y.upper == y.lower
                        yrec_convex = yrec_concave;
                    else
                        yrec_convex = (yrecL-yrecU)./(y.upper-y.lower).*mid_cv + (y.upper.*yrecU-y.lower.*yrecL)/(y.upper-y.lower);
                    end
                end
                yrec = McCormick(yrecL, yrecU, yrec_convex, yrec_concave);
            end           
            z = times(x,yrec);    
        end   
        
        %overload x.^y
        %note: y must be a positive integer
        function z = power(x,y)
            if (~isa(x,'McCormick'))
                z = McCormick(power(x,y), power(x,y), power(x,y), power(x,y));
            else       
                if mod(y,2) == 0 % when y is an even number
                    if x.lower >= 0
                        zL = x.lower.^y;
                        zU = x.upper.^y;
                        Fcv_min = x.lower;
                        Fcc_max = x.upper;
                    end
                    if x.upper <= 0
                        zL = x.upper.^y;
                        zU = x.lower.^y;
                        Fcv_min = x.upper;
                        Fcc_max = x.lower;
                    end
                    if (x.lower < 0) && (x.upper > 0)
                        zL = 0;
                        Fcv_min = 0;
                        zU = max(x.lower.^y, x.upper.^y);
                        if abs(x.lower) >= abs(x.upper)
                            Fcc_max = x.lower;
                        else
                            Fcc_max = x.upper;
                        end
                    end  
                    mid_cv = median([x.convex,x.concave,Fcv_min]);
                    z_convex = mid_cv.^y;
                    mid_cc = median([x.convex,x.concave,Fcc_max]);
                    if x.lower == x.upper
                        z_concave = z_convex;
                    else
                        z_concave = (x.lower.^y-x.upper.^y)./(x.lower-x.upper).*mid_cc + (x.lower.*x.upper.^y-x.upper.*x.lower.^y)./(x.lower-x.upper);
                    end
                else %when y is an odd number
                    zL = x.lower.^(y);
                    zU = x.upper.^(y);
                    Fcv_min = x.lower;
                    Fcc_max = x.upper;
                    mid_cv = median([x.convex,x.concave,Fcv_min]);
                    mid_cc = median([x.convex,x.concave,Fcc_max]); 
                    %convex and concave envelopes of odd degree monomials by Liberti and Pantelides (2003) are employed 
                    k = (y-1)./2;
                    a = x.lower;
                    b = x.upper;
                    r = [-0.5, -0.6058295862, -0.6703320476, -0.7145377272, -0.7470540749, -0.7721416355, -0.7921778546, -0.8086048979, -0.8223534102, -0.8340533676];
                    c = r(k).*a;
                    d = r(k).*b;
                    R = (r(k).^(2.*k+1)-1)./(r(k)-1);
                    if c < b
                        if mid_cv < c
                            z_convex = a.^(2.*k+1).*(1+R.*(mid_cv./a-1));
                        else
                            z_convex = mid_cv.^(2.*k+1);
                        end
                    else
                        z_convex = a.^(2.*k+1)+(b.^(2.*k+1)-a.^(2.*k+1))./(b-a).*(mid_cv-a);
                    end
                    if d > a
                        if mid_cc > d
                            z_concave = b.^(2.*k+1).*(1+R.*(mid_cc./b-1));
                        else
                            z_concave = mid_cc.^(2.*k+1);
                        end
                    else
                        z_concave = a.^(2.*k+1)+(b.^(2.*k+1)-a.^(2.*k+1))./(b-a).*(mid_cc-a);
                    end
                end
                z = McCormick(zL,zU,z_convex,z_concave);               
            end  
        end
        %overload x + y
        function z = plus(x,y)
             if (~isa(x,'McCormick'))
                x = McCormick(x,x,x,x);
            end
            if (~isa(y,'McCormick'))
                y = McCormick(y,y,y,y);
            end 
            z = McCormick(x.lower+y.lower, x.upper+y.upper, x.convex+y.convex, x.concave+y.concave);
        end
        %overload -x
        function z = uminus(x)
            if (~isa(x,'McCormick'))
                z = McCormick(-x, -x, -x, -x);
            else
            z = McCormick(-x.upper, -x.lower, -x.concave, -x.convex);
            end
        end
        %overload x - y
        function z = minus(x,y)
            z = plus(x,uminus(y));
        end
    end
end
      
            
            
            
            
            
      