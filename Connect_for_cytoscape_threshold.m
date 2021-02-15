%   Connect_for_cytoscape is function used to turn the network from matrix
%   into list format which can be imported to the cytoscape software to
%   draw graph networks. The output of function Connect_for_cytoscape is
%   three clumn list: first column is TFs,second column is target genes, 
%   third coulumn is regulatory strenghness.  
%   Threshold is parameter to control the scale of the network by
%   deleleting the low regulatory weightness.  
 
function [A_result]=Connect_for_cytoscape_threshold(threshold,J,name_TF,name_gene)
% threshold is decide to retain the percentage of the edges


[rows1,cols1,vals1] = find(J>=threshold);
[rows2,cols2,vals2] = find(J<=-threshold);
B=[rows1,cols1,vals1];C=[rows2,cols2,vals2];
A=[B;C];
 
A_first=A(:,2);
A_second=A(:,1);

%[n,m]=size(A);
n=size(A,1);
 
A_first_sym1=cell(1,n);A_second_sym1=cell(1,n);A_third=cell(1,n);
% AA_second_1=[];
for i=1:n
     A_first_sym1(i)=name_TF(A_first(i));
     A_second_sym1(i)=name_gene(A_second(i,1));
     A_third(i)=num2cell(J(A_second(i),A_first(i)));
      
end
A_first_sym=A_first_sym1';
A_second_sym=A_second_sym1';
A_third=A_third';
% A_third=num2cell(A_third);
A_result=[A_first_sym,A_second_sym,A_third];
end

 