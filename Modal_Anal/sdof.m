% Method From page 6-77 Allemang
 % These from vtb7_4 example                                                     
% -5.00e-03 + 1.00e+00i     5.00e-03       1.00e+00    
% -5.00e-03 - 1.00e+00i     5.00e-03       1.00e+00    
% -1.50e-02 + 1.73e+00i     8.66e-03       1.73e+00    
% -1.50e-02 - 1.73e+00i     8.66e-03       1.73e+00    
%w=f2*2*pi*i;
%w=w(110:120);
%R=-R2;%-R2(110:120);

% Above here is note for Greg

i=sqrt(-1);



w=[30;40;50]*i
R=[(7.1298e-5-3.0556e-6*i);  -1.9109e-18-1.25e-3*i; -5.5385e-5-3.0769e-6*i]





c=-w.^2.*R;
a0=w.^0.*R;
a1=w.^1.*R;
a2=R*0-1;
a=[a0 a1 a2];
%a(1:10,:)
b=a\c;
rs=roots([1 b(2) b(1)])
f=abs(rs(1))
z=real(rs(1))/abs(rs(1))
damp(rs)
b
disp('Note the wacky root (rs(1)). This is a residual, I think')
disp('Now we include the negative (j\omega) part of the axis to get the conjugates in the answer)')
pause


ll=length(w)
w2=w*0;
R3=R*0;
for i=1:ll
	R3(i)=conj(R(ll+1-i));
	w2(i)=conj(w(ll+1-i));
end
w=[w2;w];
R=[R3;R];

c=-w.^2.*R;
a0=w.^0.*R;
a1=w.^1.*R;
a2=R*0-1;
a=[a0 a1 a2];
%a(1:10,:)
b=a\c
rs=roots([1 b(2) b(1)])
f=abs(rs(1))
z=real(rs(1))/abs(rs(1))
damp(rs)

