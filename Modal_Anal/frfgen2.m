function H=frfgen2(A0,A1,B0,B1,B2,f)

om=2*pi*j*f';
H=[];
I=eye(size(A0));
for i=1:length(om)
    Hnew=-((A0+A1*om(i)+I*om(i)^2)\(B0+B1*om(i)+B2*om(i)^2)).';
   
    H=[H;Hnew];
end
