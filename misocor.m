function [x,tcm] = misocor(z,formula,tracer,purity)
%--------------------------------------------------------------------------
%Note

%This is a modified/MATLAB version of IsoCor. This script is created by
%Daniel Du from Proteomics and Metabolimics Core Facility at MD Anderson 
%Cancer Center. Anyone uses this script for publication purposes should 
%cite the original IsoCor application note, for their substantial 
%improvement on the construction of correction matrix. The original IsoCor 
%application note can be found here.

%PMID: 22419781 DOI: 10.1093/bioinformatics/bts127

%--------------------------------------------------------------------------
%Parameter Definition

%z, fraction abundance of measured ion (FAM), e.g. [0.9 0.1 0 0 0 0 0]
%formula, the chemical formula of the compound of interest, e.g. C6H12O6
%tracer, the tracer element, e.g. C
%x, mass distribution vector (MDV) with unwanted contribution removed
%cm, the correction matrix

%--------------------------------------------------------------------------
%% 0 isotope data
el_lib = {'O','H','N','C','Si','S'}; %library element
iso_dist_lib = {[0.99757	0.00038	0.00205],[0.99985	0.00015],[0.99632...	
    0.00368],[0.9893	0.0107],[0.922297	0.046832	0.030872],...
    [0.9493	0.0076	0.0429	0	0.0002]};%natural abundances for the library elements

%% 1 parse formula
parsed_formula = parse_formula(formula); %parse chemical formula into a struct
find_tracer = find(strcmp(fieldnames(parsed_formula),tracer)==1,1); %locate tracer atom in the chemical formula
if isempty(find_tracer) == 1 
    error('The tracer element is not found in the formula!\n'); 
    %report error if the tracer atom is not present in the chemical formula
else
    nt = getfield(parsed_formula,tracer); %the number of tracer atoms, and different from Nt
end
parsed_formula_q = rmfield(parsed_formula,tracer); %parsed formula for non-tracer only
nontracer = fieldnames(parsed_formula_q); %non-tracer atom list, similar to tracer

%% 2 construct total correction matrix
if iscolumn(z) == 0; z = z'; end; %make sure z is a column vector
mdv_length = length(z);
if mdv_length >= (nt+1)
    z = z(1:nt+1); %FAM beyond nt+1 components makes no contribution to MDV of the tracer, and is thus truncated
else
    z = [z;zeros(nt+1-mdv_length,1)]; %pad zeros if FAM is not long enough
end
z = z/sum(z); %normalize z
tcm = eye(nt+1); %initialize the total correction matrix

%CM of non-tracer element(s)
if isempty(nontracer) == 0
    for k = 1:length(nontracer)
        el = nontracer{k}; %the non-tracer element
        nq = getfield(parsed_formula,el); %the number of non-tracer atom, and different from Nq
        el_index = find(strcmp(el_lib,el)==1,1); %calculate correction matrix only when the element is in the library
        if isempty(el_index) == 0
            iso_dist = iso_dist_lib{el_index}; %natural abundance of the non-tracer element
            cv = 1;
            for j = 1:nq
                cv = conv(cv,iso_dist);
            end
            cm = zeros((nt+nq*(length(iso_dist)-1)+1),(nt+1)); %correction matrix
            pcv = [cv';zeros(nt,1)]; %padded correction vector
            cm(:,1) = pcv;
            for i = 2:nt+1
                cm(:,i) = cm([end 1:end-1],i-1); %construct the column vector
            end
            tcm = tcm*cm(1:nt+1,:); %truncate the rows beyond nt+1
        end
    end
end

%CM of tracer element
cm = zeros(nt+1);
iso_dist = iso_dist_lib{strcmp(el_lib,tracer)==1}; %natural abundance of the tracer element
for i = 1:nt+1
    cv = 1;
    for j = 1:nt+1-i
        cv = conv(cv,iso_dist); %use convolution to obtain a column vector of the matrix
    end
    cm(:,i) = [zeros(i-1,1);cv']; %pad the column vector with zeros
end
tcm = tcm*cm;

%CM of impurity
cm = zeros(nt+1);
for i = 1:nt+1
    cv = 1;
    for j = 1:i-1
        cv = conv(cv,[1-purity purity]); %use convolution to obtain a column vector of the matrix
    end
    cm(:,i) = [cv';zeros(nt+1-i,1)]; %pad the column vector with zeros
end
tcm = tcm*cm;

%% 3 constrained regression
function f = misocor_cost(x)
f = norm(tcm*x-z); %cost function
end
x0 = ones(nt+1,1); %initial guess
x = fmincon(@misocor_cost,x0,[],[],ones(1,nt+1),1,zeros(nt+1,1),[]); 
%optimize x, here normalization is imposed as a constraint
%Aeq = ones(nt+1,1), and beq = 1, means sum is unity
x = x(1:min(mdv_length,nt+1)); %truncate MDV to the input length if input length is shorter

end