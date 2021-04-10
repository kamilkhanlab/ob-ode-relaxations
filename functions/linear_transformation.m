%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A function that defines the linear transformation "v" in the article.    %
%                                                                         %
%Last modified by Yingkai Song 08/18/2020                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = linear_transformation(alpha,xicv,xicc)
v = (alpha+1)./2.*xicc-(alpha-1)./2.*xicv;
end