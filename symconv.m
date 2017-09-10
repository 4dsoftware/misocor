function f = symconv(a,b)
%--------------------------------------------------------------------------
%Note

%This script is created by Daniel Du from Proteomics and Metabolimics Core 
%Facility at MD Anderson Cancer Center. 

%This script calculates convolution of symbols rather than values.
%--------------------------------------------------------------------------
la = length(a);
lb = length(b);
f = sym(zeros(1,la+lb-1));
for k = 1:la+lb-1
    for j = 1:k
        if j<=la && k+1-j<=lb
            f(k) = f(k) + a(j)*b(k+1-j);
        end
    end
end