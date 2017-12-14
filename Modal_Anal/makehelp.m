cd 'Documents:MyMath:vibrate:'
pwd
filename='asd.m';

dn=0;
fidr=fopen(filename,'rt');
fidw=fopen(filename,'w');
l=fgets(fidr);
fprintf(fidw,[ l(1:length(l)) ]);
	
	
while dn~=1
	l=fgets(fidr);
	if strcmp(l(1),'%')
		fprintf(fidw,['%' l(1:length(l)) ]);
	else
		dn=1;
	end
end
dn=0;

l=fgets(fidr);
fprintf(fidw,'\n');

	
while dn~=1
	l=fgets(fidr);
	if strcmp(l(1),'%')
		fprintf(fidw,['%' l(1:length(l)) ]);
	else
		dn=1;
	end
end


	fclose(fidw);	
	fclose(fidr);	
pcode asd


