clear all
close all
clc

dbstop if error
profile on

a=testSkel();

profile viewer

% remove incorrectly identified nonrunnable lines (end statements after
% break, continue, return and error statements)
nonrunnable={
    'Sskel',7
    'Sskel>sortInto',1
    'checkManifold',3
    'inset',2
    'closure',2
    'deSliver',1
    'polyCompare',2
    'vertex>vertex.vertex',2
    };

stats=profile('INFO');
cov=internal.matlab.codetools.reports.buildCoverageInfo(cd,'temp');
k=0;
for i=1:length(cov)
    if isempty(cov(i).funlist)
        continue
    end
    for j=1:length(cov(i).funlist)
        if strcmp(cov(i).funlist(j).name,'testSkel')
            continue
        end
        k=k+1;
        coverage{k,2}=cov(i).funlist(j).name;
        linestotal(k)=cov(i).funlist(j).runnablelines;
        linesrun(k)=round(linestotal(k)*cov(i).funlist(j).coverage/100);
        ind=strcmp(coverage{k,2},nonrunnable(:,1));
        m=find(ind);
        if ~isempty(m)
            linestotal(k)=linestotal(k)-nonrunnable{m,2};
        end
        coverage{k,1}=linesrun(k)/linestotal(k)*100;
    end
end
coverage
totalLines=sum(linestotal)
totalCoverage=sum(linesrun)/totalLines*100