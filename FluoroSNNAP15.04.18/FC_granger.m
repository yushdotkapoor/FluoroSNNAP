function [sig,pval,F] = FC_granger(dF_cell)
load('params.mat');
X = dF_cell;
momax = params.FC.GC.morder; % Maximum model order
icregmode = 'LWR';
regmode = 'OLS';
mhtc = 'FDR';
% for some reason, the value of alpha is morder???
alpha = params.FC.GC.morder;
% disp(X)
% disp(momax)
% disp(icregmode)
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
morder = 1;
[A,SIG] = tsdata_to_var(X,morder,regmode);
assert(~isbad(A),'VAR estimation failed');

%%
if(params.FC.GC.iter==-1)
[G,info] = var_to_autocov(A,SIG,[]);
else
    [G,info] = var_to_autocov(A,SIG,params.FC.GC.iter);
end
var_info(info,true); 
F = autocov_to_pwcgc(G);

%%
assert(~isbad(F,false),'GC calculation failed');
pval = mvgc_pval(F,morder,size(dF_cell,2),1,1,1,size(dF_cell,1)-2,'F');
sig  = significance(pval,alpha,mhtc);
sig(isnan(sig))=0;
sig = sig';
pval = pval';
F = F';
