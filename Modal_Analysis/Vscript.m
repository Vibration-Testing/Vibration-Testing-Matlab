function V=Vscript(H,omega,plist)

% V=Vscript(H, omega, plist)
%
% Returns the variable script V defined by Slater in "Vibration Testing". 
%
% H is the matrix of FRFs, in blocks of p outputs by o inputs, stacked
%   vertically by increasing frequency for each block.
% omega is the list of frequencies (real) in radians per second
% plist are the powers of the blocks output. 
%   0 would return H with the blocks each transposed
%   2 would return H with each transposed block multiplied by jw^2
%   [0 2] would return the two block columns of the preceeding examples
%   combined in one matrix.
%
% 

global freqdebug
if freqdebug==1, disp('Making Vscript'),end
p=size(H,1)/length(omega);
o=size(H,2);
if floor(p)~=p
    disp(['On cannot have ' num2str(size(H,1)/length(omega)) ' outputs'])
return
else
    H2=zeros(o*length(omega),p);
    if freqdebug==1, disp('Regigering H'),end
    for ii=1:length(omega)
        H2(((ii-1)*o+1):ii*o,:)=transpose(H(((ii-1)*p+1):ii*p,:));
    end
    if freqdebug==1, disp('Clearing H'),end
    
    clear('H')
    if freqdebug==1, disp('Making empty V'),end

    V=zeros(size(H2,1),size(H2,2)*length(plist));

    % Diagonal matrix of proper size to do multiplications in one fell swoop
    % for each block column.
    if freqdebug==1, disp('Making ODV'),end 
    ODV=sparse(zeros(size(H2,1),1));
    no=length(omega);
    onev=ones(o,1);
    if freqdebug==1, disp('Putting values into ODV'),tic,end
    
    for ii=1:no
        ODV(((ii-1)*o+1):(ii*o),1)=omega(ii)*sqrt(-1)*onev;
    end
    if freqdebug==1, disp('Putting values into ODV'),toc,end
    
    

    if freqdebug==1, disp('Make OD from sparse ODV'),end
    OD=(diag(ODV));
    if freqdebug==1, disp('Finish making V'),end
    for ii=1:length(plist)
        V(:,((ii-1)*p+1):(ii*p))=(OD^plist(ii))*H2;
    end
    
    
end