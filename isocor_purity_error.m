%--------------------------------------------------------------------------
%Note

%This script is created by Daniel Du from Proteomics and Metabolimics Core 
%Facility at MD Anderson Cancer Center. 

%This script compares the IsoCor approach to construct the a single matrix 
%for both natural abundance and isotopic impurity and the standard approach
%to construct matrices for natural abundance and isotopic impurity 
%respectively. The software IsoCor can be found here.

%PMID: 22419781 DOI: 10.1093/bioinformatics/bts127

%--------------------------------------------------------------------------
%Parameter Definition

%tracer:C
%natural abundance: 1-a, a
%purity: b is purity of 13C in the nutrient 

%--------------------------------------------------------------------------

clear
clc
nc = 5; %number of carbon atoms
syms a b %initialize symbols

%% 1 constructing matrices for natural abundance and isotopic impurity respectively (standard)
%correction matrix for natural abundance
cm = sym(zeros(nc+1)); %initialize the matrix
for i = 1:nc+1
    tmp = 1;
    for j = 1:nc+1-i
        tmp = symconv(tmp,sym([1-a a])); %use convolution to obtain a column vector of the matrix
    end
    cm(:,i) = [zeros(i-1,1);transpose(tmp)]; %pad the column vector with zeros
end

%correction matrix for isotopic impurity
pm = sym(zeros(nc+1));
for i = 1:nc+1
    tmp = 1;
    for j = 1:i-1
        tmp = symconv(tmp,sym([1-b b])); %use convolution to obtain a column vector of the matrix
    end
    pm(:,i) = [transpose(tmp);zeros(nc+1-i,1)]; %pad the column vector with zeros
end

%matrix multiplication to obtain the total correction matrix
mat1 = pm*cm;

%% 2 constructing a single matrix for both natural abundance and isotopic impurity (IsoCor)
mat2 = sym(zeros(nc+1));
for i = 1:nc+1
    tmp = 1;
    for j = 1:nc+1-i
        tmp = symconv(tmp,sym([1-a a])); %use convolution to obtain a column vector of the correction matrix for natural abundance
    end
    for j = 1:i-1
        tmp = symconv(tmp,sym([1-b b])); %use convolution to obtain a column vector of the correction matrix for isotopic impurity
    end
    mat2(:,i) = transpose(tmp); 
end

%% 3 Compare
%expand the expression to remove parentheses
expand(mat1 - mat2)

%conclusion: the IsoCor approach, despite the concise format, is 
%oversimplified and not equivalent to the standard approach. 
