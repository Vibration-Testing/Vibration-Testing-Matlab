function U=Uscript(o,omega,plist)

% U=Uscript(o, omega, plist)
% o is number of inputs
% omega is frequency vector in radians per second
% plist is a vector of powers of omega in adjacent columns 
global freqdebug
if freqdebug==1, disp('Making Vscript'),end
U=zeros(o*length(omega),o*length(plist));
Ip=eye(o);
no=length(omega);
onev=ones(o,1);
for ii=1:length(omega)
    for jj=1:length(plist)
        U(((ii-1)*o+1):ii*o,((jj-1)*o+1):jj*o)=eye(o)*(sqrt(-1)*omega(ii))^plist(jj);
    end
end
