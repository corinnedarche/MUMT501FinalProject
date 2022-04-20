function [dFunNorm] = postprocess(dFun, gamma)
% FUNCTION: postprocess.m
% Created by: Corinne Darche
% This function applies the postprocessing step for onsetImplementation.m
% to each detection function

dFunNorm = dFun - mean(dFun(:,:));
dFunNorm = dFunNorm./std(dFun(:,:));

b = [1 -1*gamma];
a = [1 -1*gamma];

dFunNorm = filter(b,a,dFunNorm);

end

