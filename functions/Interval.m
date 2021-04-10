%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Natural interval extension via operator overloading                      % 
%                                                                         %
%Overloaded operators: +, -, ./, .*, exp, log, x.^n                       % 
%                                                                         %
%Last modified by Yingkai Song, 08/17/2020                                %
%                                                                         %
%Reference                                                                %
%          Moore, R.E.: Methods and Applications of Interval Analysis.    % 
%          SIAM, Philadelphia (1979)                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Interval
    properties
        lower %lower bound
        upper %upper bound
    end
    methods
        function z = Interval(lower, upper)
            z.lower = lower;
            z.upper = upper;
        end
        %overload x + y
        function z = plus(x, y)
            %turn non-Interval object to Interval object
            if ~isa(x, 'Interval')
                x = Interval(x, x);
            end
            if ~isa(y, 'Interval')
                y = Interval(y, y);
            end
            z = Interval(x.lower+y.lower, x.upper+y.upper);
        end
        %overload -x
        function z = uminus(x)
            if ~isa(x, 'Interval')
                z = Interval(-x, -x);
            else
                z = Interval(-x.upper, -x.lower);
            end 
        end
        %overload x - y
        function z = minus(x, y)
            z = plus(x, uminus(y));
        end
        %overload x .* y
        function z = times(x, y)
            if ~isa(x, 'Interval') 
                x = Interval(x, x);
            end
            if ~isa(y, 'Interval')
                y = Interval(y, y);
            end
            t = [x.lower.*y.lower, x.lower.*y.upper, x.upper.*y.lower, x.upper.*y.upper];
            z = Interval(min(t), max(t));
        end
        
        %overload x ./ y
        %0 cannot be in [y.lower,y.upper]
        function z = rdivide(x, y)
            %this is overloaded by x.*(1./y)
            if ~isa(y, 'Interval')
                yrec = Interval(1./y, 1./y);
            else
                yrec = Interval(1./y.upper, 1./y.lower);
            end       
            z = times(x, yrec);    
        end     
        
        %overload exp(x)
        function z = exp(x)
            if ~isa(x, 'Interval')
                z = Interval(exp(x), exp(x));
            else
                z = Interval(exp(x.lower), exp(x.upper));
            end  
        end
        
        %overload log(x)
         function z = log(x)
            if ~isa(x, 'Interval')
                z = Interval(log(x), log(x));
            else
                z = Interval(log(x.lower), log(x.upper));
            end
         end
        
        %overload x.^y
        %note: y must be a positive integer
        function z = power(x,y)
            if ~isa(x, 'Interval')
                z = Interval(power(x,y), power(x,y));
            else
                if mod(y,2) %when y is an odd number 
                   zL = x.lower.^y;
                   zU = x.upper.^y;
                else %when y is an even number
                    zU = max(x.lower.^y, x.upper.^y);
                    if (x.lower >= 0) || (x.upper <= 0)
                        zL = min(x.lower.^y, x.upper.^y);  
                    else
                        zL = 0;
                    end
                end
                z = Interval(zL, zU);   
            end
        end 
        
    end
end
              
            
            
            