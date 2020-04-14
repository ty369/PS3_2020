% ChemE 7770 PS3 #2 Check Balance
%bound=xlsread('SMatrix.xlsx','S','H2:V19');
%Z=transpose(Element)*bound

S=xlsread('SMatrix.xlsx','S','B2:V19');
Element=xlsread('SMatrix.xlsx','Elements');
Rxn=xlsread('SMatrix.xlsx','S','B2:G19'); %checek balance on the 6 rxn
X=transpose(Element)*Rxn; %Inside the box
Y=transpose(Element)*S; % Outside the box
display(X);
display(Y);
% Internal reactions are balanced. 

