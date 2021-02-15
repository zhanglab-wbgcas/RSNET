%   WeakVarianceDelete is function to delete the genes and related gene
%   expression profils with low variance.Threshold is parameter used to control the
%   strenthness of variance;Textdata is gene list; Data is gene expression profile
%   Parameter can be set according to the requirment and dataset,for
%   example value between 0.05 and 0.5 
function [textdata,data]=WeakVarianceDelete(textdata,data,threshold)
vardata=var(data,0,2); % computer the variance of row
%[m,n1]=size(vardata);
m=size(vardata,1);
k=0;
for i=1:m
    if vardata(i,1)<threshold
       k=k+1;
       data(1-k+i,:)=[];
       textdata(1-k+i,:)=[];
    end
end
end

