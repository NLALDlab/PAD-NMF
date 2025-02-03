close all;
clear;
clc;

%% ------------------------------------------------------------------------------------------------------------------------------------%%
%Configuration
BlocksNum = 3;
BlocksSize = 8*ones(1,BlocksNum);
CSBlockSize = [0,cumsum(BlocksSize)];

%Assemble threshold vector
load('TestingOutput/TestingSetPianoOutput_Dataset3_PADNMF_Training.mat');

InstancesNum = size(TestingSetOutput,2);

GeoMeans = zeros(BlocksNum,InstancesNum);

for i = 1:InstancesNum
    HCurr = TestingSetOutput{4,i}{end};
    for j = 1:BlocksNum
        HBlock = HCurr(1+CSBlockSize(j):CSBlockSize(j)+BlocksSize(j),:);
        
        Temp = [];
        for k = -(BlocksSize(j)-1):(BlocksSize(j)-1)
            Temp = [Temp, GeoMeanVec(diag(HBlock,k))]; 
        end

        GeoMeans(j,i) = max(Temp);
    end
end

Thresh = max(GeoMeans(:,[4,8,12]),[],2);


%Output to show
load('TestingOutput/TestingSetPianoOutput_Dataset3_PADNMF_Testing.mat');

InstancesNum = size(TestingSetOutput,2);

GeoMeans = zeros(BlocksNum,InstancesNum);

for i = 1:InstancesNum
    HCurr = TestingSetOutput{4,i}{end};
    for j = 1:BlocksNum
        HBlock = HCurr(1+CSBlockSize(j):CSBlockSize(j)+BlocksSize(j),:);
        
        Temp = [];
        for k = -(BlocksSize(j)-1):(BlocksSize(j)-1)
            Temp = [Temp, GeoMeanVec(diag(HBlock,k))]; 
        end

        GeoMeans(j,i) = max(Temp);
    end
end

%Show geometric means
symb = {'*','o','square'};
colors = {'r','m','b'};

ColLow = 0.3;
ColHigh = 0.8;

%Labels for first 30 outputs 
NoteTable = [1,0,0,0,2,0,3,0,0,0,1,0,0,2,0,3,0,0,0,0,1,0,0,2,0,0,3,0,0,0];
NoiseTable = ~NoteTable;

figure(1);

for i = 1:30
    ColList = linspace(ColLow,ColHigh,BlocksNum);
    for j = 1:BlocksNum
        if j == NoteTable(i)
            COL = [0,0,0];
            MS = 9;
            LW = 2;
        else
            COL = ColList(1)*ones(1,BlocksNum);
            MS = 7;
            LW = 1;
        end

        plot(i,10*log10(GeoMeans(j,i)),symb{j},'Color',COL,'MarkerSize',MS,'LineWidth',LW);
        hold on;
        ColList = ColList(2:end);
    end
end

%Plot thresholds
ColList = linspace(ColLow,ColHigh,BlocksNum);
for i = 1:BlocksNum
    yline(10*log10(Thresh(i)),'-.','Color',ColList(i)*ones(1,BlocksNum),'LineWidth',1);
end

xlabel('Mixture number','FontSize',18,'Position',[15,-166]);
ylabel('Geometric mean (dB)','FontSize',18,'Position',[-0.8,-90]);
grid on;
grid minor;


function GM = GeoMeanVec(v)
    split = floor(length(v)/3);
    v1 = v(1:split);
    v2 = v(split+1:2*split);
    v3 = v(2*split+1:end);
    GM = ( prod(v1)^(1/length(v)) )*( prod(v2)^(1/length(v)) )*( prod(v3)^(1/length(v)) ); %to avoid underflow
end