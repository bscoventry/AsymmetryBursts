close all;clear all;
load('AsymBurstData.mat')
numSubjects = 24;
numConditions = 4;
stimulusCondA = cell(1,1);
stimulusCondV = cell(1,1);
subject = cell(1,1);
totalCountsData = [];
totCountsSubject = cell(1,4);
totSubjectClass = cell(1,4);
counter = 1;
for ck = 1:numConditions
    curCount = totData{2,ck};
    if ck == 1
        ACond = 'A2';
        VCond = 'V2';
    elseif ck == 2
        ACond = 'A8';
        VCond = 'V2';
    elseif ck == 3
        ACond = 'A2';
        VCond = 'V8';
    elseif ck == 4
        ACond = 'A8';
        VCond = 'V8';
    end
    for bc = 1:numSubjects
        curData = curCount{1,bc};
        dataLen = length(curData);
        for jk = 1:dataLen
            totalCountsData = [totalCountsData curData(jk)];
            stimulusCondA{counter} = ACond;
            stimulusCondV{counter} = VCond;
            subject{counter} = num2str(bc);
            counter = counter + 1;
        end
    end
end
hTotCounts = kstest(totalCountsData);
p = anovan(totalCountsData,{stimulusCondA,stimulusCondV},'model','interaction','varnames',{'StimulusCondA','StimulusCondV'});
for ck = 1:numConditions
    curCount = totData{2,ck};
    totSubData = [];
    subjectClass = [];
    for bc = 1:numSubjects
        curData = curCount{1,bc};
        dataLen = length(curData);
        for jk = 1:dataLen
            totSubData = [totSubData curData(jk)];
            subjectClass = [subjectClass ck];
        end
    end
    totCountsSubject{ck} = totSubData;
    totSubjectClass{ck} = subjectClass;
end
ConditionList = {'V2A2','V2A8','V8A2','V8A8'};
figure(2)
C = [totCountsSubject{1} totCountsSubject{2} totCountsSubject{3} totCountsSubject{4}];
L = [totSubjectClass{1} totSubjectClass{2} totSubjectClass{3} totSubjectClass{4}];
boxplot(C,L,'Notch','on','Labels',{'V2A2','V2A8','V8A2','V8A8'})
%Low Count
stimulusCondA = cell(1,1);
stimulusCondV = cell(1,1);
totalCountsLowData = [];
counter = 1;
for ck = 1:numConditions
    curCount = totData{3,ck};
    if ck == 1
        ACond = 'A2';
        VCond = 'V2';
    elseif ck == 2
        ACond = 'A8';
        VCond = 'V2';
    elseif ck == 3
        ACond = 'A2';
        VCond = 'V8';
    elseif ck == 4
        ACond = 'A8';
        VCond = 'V8';
    end
    for bc = 1:numSubjects
        curData = curCount{bc,1};
        dataLen = length(curData);
        for jk = 1:dataLen
            totalCountsLowData = [totalCountsLowData curData(jk)];
            stimulusCondA{counter} = ACond;
            stimulusCondV{counter} = VCond;
            counter = counter + 1;
        end
    end
end
hTotLowCounts = kstest(totalCountsLowData);
pLow = anovan(totalCountsLowData,{stimulusCondA,stimulusCondV},'model','interaction','varnames',{'StimulusCondA','StimulusCondV'});
figure()
for ck = 1:numConditions
    curCount = totData{3,ck};
    totSubLowData = [];
    subjectClass = [];
    for bc = 1:numSubjects
        curData = curCount{bc,1};
        dataLen = length(curData);
        for jk = 1:dataLen
            totSubLowData = [totSubLowData curData(jk)];
            subjectClass = [subjectClass ck];
        end
    end
    totCountsSubjectLow{ck} = totSubLowData;
    totSubjectClassLow{ck} = subjectClass;
end
CLow = [totCountsSubjectLow{1} totCountsSubjectLow{2} totCountsSubjectLow{3} totCountsSubjectLow{4}];
LLow = [totSubjectClassLow{1} totSubjectClassLow{2} totSubjectClassLow{3} totSubjectClassLow{4}];
boxplot(CLow,LLow,'Notch','on','Labels',{'V2A2','V2A8','V8A2','V8A8'})

%High Count
stimulusCondA = cell(1,1);
stimulusCondV = cell(1,1);
totalCountsHighData = [];
counter = 1;
for ck = 1:numConditions
    curCount = totData{4,ck};
    if ck == 1
        ACond = 'A2';
        VCond = 'V2';
    elseif ck == 2
        ACond = 'A8';
        VCond = 'V2';
    elseif ck == 3
        ACond = 'A2';
        VCond = 'V8';
    elseif ck == 4
        ACond = 'A8';
        VCond = 'V8';
    end
    for bc = 1:numSubjects
        curData = curCount{bc,1};
        dataLen = length(curData);
        for jk = 1:dataLen
            totalCountsHighData = [totalCountsHighData curData(jk)];
            stimulusCondA{counter} = ACond;
            stimulusCondV{counter} = VCond;
            counter = counter + 1;
        end
    end
end
hTotHighCounts = kstest(totalCountsHighData);
pHigh = anovan(totalCountsHighData,{stimulusCondA,stimulusCondV},'model','interaction','varnames',{'StimulusCondA','StimulusCondV'});

figure()
for ck = 1:numConditions
    curCount = totData{4,ck};
    totSubHighData = [];
    subjectClass = [];
    for bc = 1:numSubjects
        curData = curCount{bc,1};
        dataLen = length(curData);
        for jk = 1:dataLen
            totSubHighData = [totSubHighData curData(jk)];
            subjectClass = [subjectClass ck];
        end
    end
    totCountsSubjectHigh{ck} = totSubHighData;
    totSubjectClassHigh{ck} = subjectClass;
end
CLow = [totCountsSubjectHigh{1} totCountsSubjectHigh{2} totCountsSubjectHigh{3} totCountsSubjectHigh{4}];
LLow = [totSubjectClassHigh{1} totSubjectClassHigh{2} totSubjectClassHigh{3} totSubjectClassHigh{4}];
boxplot(CLow,LLow,'Notch','on','Labels',{'V2A2','V2A8','V8A2','V8A8'})
%Calculate High Powers

stimulusCondA = cell(1,1);
stimulusCondV = cell(1,1);
PowerHighData = [];

counter = 1;
for ck = 1:numConditions
    curCount = totData{5,ck};
    if ck == 1
        ACond = 'A2';
        VCond = 'V2';
    elseif ck == 2
        ACond = 'A8';
        VCond = 'V2';
    elseif ck == 3
        ACond = 'A2';
        VCond = 'V8';
    elseif ck == 4
        ACond = 'A8';
        VCond = 'V8';
    end
    for bc = 1:numSubjects
        curData = curCount{1,bc};
        dataLen = length(curData);
        for jk = 1:dataLen
            curPow = curData{1,jk};
            for tk = 1:length(curPow)
                curPulse = curPow{1,tk};
                if ~isinf(mean(curPulse))
                    PowerHighData = [PowerHighData mean(curPulse)];
                    stimulusCondA{counter} = ACond;
                    stimulusCondV{counter} = VCond;
                    counter = counter + 1;
                end
            end
        end
    end
end
hPowHighCounts = kstest(PowerHighData);
pPowHigh = anovan(PowerHighData,{stimulusCondA,stimulusCondV},'model','full','varnames',{'StimulusCondA','StimulusCondV'});

%Calculate Low Powers

stimulusCondA = cell(1,1);
stimulusCondV = cell(1,1);
PowerLowData = [];

counter = 1;
for ck = 1:numConditions
    curCount = totData{6,ck};
    if ck == 1
        ACond = 'A2';
        VCond = 'V2';
    elseif ck == 2
        ACond = 'A8';
        VCond = 'V2';
    elseif ck == 3
        ACond = 'A2';
        VCond = 'V8';
    elseif ck == 4
        ACond = 'A8';
        VCond = 'V8';
    end
    for bc = 1:numSubjects
        curData = curCount{1,bc};
        dataLen = length(curData);
        for jk = 1:dataLen
            curPow = curData{1,jk};
            for tk = 1:length(curPow)
                curPulse = curPow{1,tk};
                if ~isinf(mean(curPulse))
                    PowerLowData = [PowerLowData mean(curPulse)];
                    stimulusCondA{counter} = ACond;
                    stimulusCondV{counter} = VCond;
                    counter = counter + 1;
                end
            end
        end
    end
end
hPowLowCounts = kstest(PowerLowData);
pPowerLow = anovan(PowerLowData,{stimulusCondA,stimulusCondV},'model','full','varnames',{'StimulusCondA','StimulusCondV'});