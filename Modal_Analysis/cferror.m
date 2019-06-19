function cferror=cferror(x)
global w2 R2 XoF

lx=length(x);

XoF=x(lx-2)+x(lx-1)*i*w2-x(lx)*w2.^2;
for j=1:(lx/3)-1
XoF=XoF+x(3*j-2)./(-w2.^2+2*x(3*j-1)*w2*i*x(3*j)+x(3*j)^2);%+x(4)+x(5)*i*w2-x(6)*w2.^2;
end

cferror=norm(XoF-R2);
%pause
