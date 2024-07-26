% readSNODARfromASFS
D=dir('*.dat');
dd=[];
depth=[];
for n=1:length(D)
    D2=readtable(D(n).name);
    dd=[dd;D2.TIMESTAMP];
    depth=[depth;nanmedian(D2.snowd_current_Avg)*ones(size(D2.snowd_current_Avg))];
end
