%% Check OpenSees model tcl file for unconnected nodes. 
%Two nodes elements and shellMITC4 elements are considered.
%Unconnected nodes will be listed in ErrorNodes.
%WarningNodes are only connected in 1 way, need further checking. 
%Besides, remember to check 'fix' and 'constraint' DoFs. 

%By Yunxiang Ma

%%
clear all
clc
%% read file
ModelVer='PTPlus';
SeisNum=1;
Record={'LPSLE1' 0.17;'LPSLE2' 0.19;'NSLE' 0.19;'SHSLE' 0.13;'NDBE' 0.54;'NORI1' 0.55;'IVSLE' 0.14;'NORI2' 0.55;
    'LPDBE' 0.54;'SHDBE' 0.48;'LPMCE' 0.66;'NMCE' 0.76;'SHMCE' 0.68;'N12MCE' 0.89};
RecordName=Record{SeisNum,1};


fidM=fopen(['Platform6sps' RecordName ModelVer '.tcl'],'r');
Content=textscan(fidM,'%s');
Content=convertCharsToStrings(Content{1,1});

%% rearrange file data
Length=length(Content);
NodeList=[];
ElementList=[];
FixList=[];
ConstraintList=[];

for i=1:Length
    text=Content{i};
    L=strlength(text);
    if L==4
    if text=='node'
        NodeList=[NodeList; str2num(Content{i+1})];
    end
    end
    
    if L==7
    if text=='element' 
            if strlength(Content{i+1})==10
                if Content{i+1}=='shellMITC4'
                    ShellNodes=[str2num(Content{i+3}); str2num(Content{i+4}); str2num(Content{i+5}); str2num(Content{i+6})];
                    ElementList=[ElementList;ShellNodes];
                else
                    ElementList=[ElementList;str2num(Content{i+3});str2num(Content{i+4})];
                end
            else
                ElementList=[ElementList;str2num(Content{i+3});str2num(Content{i+4})];
            end
    end
    end
    
    if L==3
    if text=='fix'
        FixList=[FixList;str2num(Content{i+1})];
    end
    end
    
    if L==8
    if text=='equalDOF'
        ConsNodes=[str2num(Content{i+1}); str2num(Content{i+2})];
        ConstraintList=[ConstraintList;ConsNodes];
    end
    end
end

%% Check node occurance in element, BC, and constraints.
Nnode=length(NodeList);
CheckMat=zeros(Nnode,5);
CheckMat(:,1)=NodeList;
for i=1:Nnode
    Node=NodeList(i);
    CheckEle=find(ElementList==Node);
    EleOccur=length(CheckEle);
    Checkfix=find(FixList==Node);
    FixOccur=length(Checkfix);
    CheckCon=find(ConstraintList==Node);
    ConOccur=length(CheckCon);
    TotOccur=EleOccur+FixOccur+ConOccur;
    CheckMat(i,2:5)=[EleOccur FixOccur ConOccur TotOccur];
end

ErrorNodeNum=find(CheckMat(:,5)==0);
ErrorNodes=CheckMat(ErrorNodeNum,1);
WarningNodeNum=find(CheckMat(:,5)==1);
WarningNodes=CheckMat(WarningNodeNum,1);
WarningNodesStatus=CheckMat(WarningNodeNum,:);