% An example for using RSNET to infer GRN: inferring GRN for apple fruit develoopment
% from gene expression data. 
% Version data: Feb.,2021
clear;clc;
%% Dataset input 
filename  = 'apple_fruit_development_data.txt';
data = importdata(filename);
data_gene_expression = data.data; % size(data_gene_expression);
data_gene_name =data.textdata; % size(data_gene_name);
sample_name = data_gene_name(1);
data_gene_name(1)= []; %data_gene_name(1:5); % size(data_gene_name)
 
%% Filter the genes with low varaince.
threshold = 5;
[data_gene_name_diff,data_gene_expression_diff] = WeakVarianceDelete(data_gene_name,data_gene_expression,threshold);
% size(data_gene_expression_diff)
% size(data_gene_name_diff)
Y =log2(data_gene_expression_diff); 
fprintf('Data prepared for %d genes! \n',size(Y,1));
%% Run RSNET method on the data
lamda =  1; 
alpha = 0.1; % parameter for correlation
gama = 0.5; % parameter for prior information
beta = 0.1; % parameter for deleting the noise
t = 0.5; %  Parameter for the interation of MI and RO;  t:[0,1]


J_na = zeros(size(Y,1),size(Y,1)); J_s=J_na;

% n_gene = size(Y,1)
n_gene = 10
% for i=1:size(Y,1)
for i=1:n_gene  % Chose few genes for running time
if mod(i,50)==0
    fprintf('Network inferring for gene: i=%d of %d.\n',i,size(Y,1));
end
y = Y(i,:);    
X = [Y(1:i-1,:);Y(i+1:size(Y,1),:)];
[net,net_value]=RSNET(y',X',lamda,alpha,gama, beta,t) ; 
J_s(i,1:i-1) = net(1:i-1); J_s(i,i+1:size(Y,1))=net(i:end);
J_na(i,1:i-1) = net_value(1:i-1); J_na(i,i+1:size(Y,1))=net_value(i:end);
end
 
 Gval=J_na; Gval=abs(Gval); G = Gval;
 q=0.5;  G(G<q) = 0;  
 for i=1:size(G,1)
     for j=1:size(G,2)
         if G(i,j)>=G(j,i)
             G(j,i)=0;
         end
     end
 end
 
% sum(G>0,1);sum(G>0,2); sum(sum(G>0));
% index=find(sum(G>0,2)>100);
% gene_list = data_gene_name_diff;
% gene_list(index) 

%% output a network for cytoscape
% turn the network from matrix to column with TF|gene|edge(0 or 1) 
gene_list = data_gene_name_diff;
threshold = 0.5;
[testfile]=Connect_for_cytoscape_threshold(threshold,G,gene_list,gene_list) ;
network_size=size(testfile,1)  
fprintf('NOTICE:\nThe Size of the Inferred Network is %d.\n',network_size);
 
xlswrite('result_network',testfile); 
