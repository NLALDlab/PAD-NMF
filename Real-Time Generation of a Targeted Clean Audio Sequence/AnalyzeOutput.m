close all;
clear;
clc;

%% ----------------------------------------------CHOOSE OUTPUT---------------------------------------- %%
Config = 2; %1: ARM constraints
            %2: No ARM constraints

if Config == 1
    Append = '_ARM';
    BlocksNum = 3;
    BlocksSize = 7*ones(1,BlocksNum);
end

if Config == 2
    Append = '_NotARM';
    BlocksNum = 3;
    BlocksSize = 8*ones(1,BlocksNum);
end

Dir = 'TestingOutput/';

%% -----------------------------------------------PARAMETERS--------------------------------------- %%
%Linear model parameters
MaxItTrain = 20000;
MaxItTest = 50000;

%Choose different values for p (0.25, 0.5, 1)
pGamma = 0.25;
pPoiss = 0.25;

switch pGamma
    case 1
        AlphaGamma = 1e-3;
        ExpGamma = 0;
    case 0.5
        AlphaGamma = 1e-1;
        ExpGamma = -0.5;
    case 0.25
        AlphaGamma = 1e-1;
        ExpGamma = -0.75;
end

switch pPoiss
    case 1
        AlphaPoiss = 1e0;
        ExpPoiss = 0;
    case 0.5
        AlphaPoiss = 2.4e-11;
        ExpPoiss = -0.5;
    case 0.25
        AlphaPoiss = 1e-2;
        ExpPoiss = -0.75;
end

%% ----------------------------CONSTRUCT LINEAR MODEL ON TRAIN DATA----------------------------- %%
load([Dir, 'TestingSetPianoOutput_Training',Append,'.mat']);
InstancesNum = size(TestingSetOutput,2);

%Compute training active indices
GeoMeans = ComputeGeoMeans(TestingSetOutput,InstancesNum,BlocksSize,BlocksNum);

%Show geometric means
symb = {'*','o','square'};
colors = {'r','m','b'};

figure(1);
title('Train geometric means');
xlabel('Instances');
ylabel('GMs');

for i = 1:InstancesNum
    for j = 1:BlocksNum
        plot(i,10*log10(GeoMeans(j,i)),[symb{j},colors{j}]);
        hold on;
    end
end

grid minor;
grid on;

%Assemble labels from training set
CGamma = [];
CPoiss = [];

for i = 1:3
    CGamma = [CGamma,eye(BlocksNum),1e-2*ones(BlocksNum,1)];
    CPoiss = [CPoiss,eye(BlocksNum),1e-8*ones(BlocksNum,1)];
end

CGamma = max(2^-52,CGamma);
CPoiss = max(2^-52,CPoiss);

%Assemble transformed active indices from training set
BGeo = GeoMeans;
BGeoLog = log10(1e17*BGeo);

BGamma = BGeo;
BPoiss = BGeoLog;

%Find linear model
[AGamma,APoiss] = FindLinearModel(BGamma,CGamma,BPoiss,CPoiss,BlocksNum,MaxItTrain);

%Compute thresholds from latent indices on training set
[xGamma,xPoiss] = InvertLinearModel(AGamma,BGamma,AlphaGamma,ExpGamma,APoiss,BPoiss,AlphaPoiss,ExpPoiss,BlocksNum,InstancesNum,MaxItTest);

ThreshGamma = 5*max(xGamma(:,[4,8,12]),[],2);
ThreshPoiss = 1.5*max(xPoiss(:,[4,8,12]),[],2);

%% ------------------------------------------APPLY LINEAR MODEL ON TEST DATA------------------------------------------------- %%
load([Dir, 'TestingSetPianoOutput_Testing',Append,'.mat']);
InstancesNum = size(TestingSetOutput,2);

%Compute testing active indices
GeoMeans = ComputeGeoMeans(TestingSetOutput,InstancesNum,BlocksSize,BlocksNum);

%Assemble transformed active indices from testing set
BGeo = GeoMeans;
BGeoLog = log10(1e17*BGeo);
BGamma = BGeo;
BPoiss = BGeoLog;

%Compute latent indices on testing set
[xGamma,xPoiss] = InvertLinearModel(AGamma,BGamma,AlphaGamma,ExpGamma,APoiss,BPoiss,AlphaPoiss,ExpPoiss,BlocksNum,InstancesNum,MaxItTest);

xGamma
xPoiss

%Slightly better thresholds 
%ThreshGamma = [0.1;0.1;0.1];
%ThreshPoiss = [0.6;0.6;0.47];

%Apply thresholds to output
xGamma(xGamma < ThreshGamma) = 0;
xPoiss(xPoiss < ThreshPoiss) = 0;

xGamma
xPoiss

%Piano notes are (roughly) contained in columns of the output whose index
%is congruent to 1 (note 1, component 1), 4 (note 2, component 2), 7 (note 3, component 3) mod 10. 
%Nonzero coefficients belonging to different columns than those, are 
%(with some exceptions) false positives. Conversely, zero columns whose index 
%is congruent to 1, 4, 7 mod 10 represent false negatives. 

%------------------------------------------------------------------------------------------------------------------------------------%
function GM = GeoMeanVec(v)
    split = floor(length(v)/3);
    v1 = v(1:split);
    v2 = v(split+1:2*split);
    v3 = v(2*split+1:end);
    GM = ( prod(v1)^(1/length(v)) )*( prod(v2)^(1/length(v)) )*( prod(v3)^(1/length(v)) ); %to avoid underflow
end

function [xGamma,xPoiss] = InvertLinearModel(AGamma,BGamma,AlphaGamma,ExpGamma,APoiss,BPoiss,AlphaPoiss,ExpPoiss,BlocksNum,InstancesNum,MaxIt)
    %Initialization
    xGamma = 1e-0*ones(BlocksNum,InstancesNum);
    xPoiss = 1e-0*ones(BlocksNum,InstancesNum);
    
    for i = 1:MaxIt
        xGamma = max(2^-52, xGamma.*( ( AGamma'*(((AGamma*xGamma).^(-2)).*BGamma) )./( AGamma'*((AGamma*xGamma).^(-1)) + AlphaGamma*xGamma.^ExpGamma ) ).^(0.5) );
        xPoiss = max(2^-52, xPoiss.*( ( APoiss'*(((APoiss*xPoiss).^(-1)).*BPoiss) )./( repmat(sum(APoiss)',1,size(xGamma,2)) + AlphaPoiss*xPoiss.^ExpPoiss ) ) );
    end
end

function [AGamma,APoiss] = FindLinearModel(BGamma,CGamma,BPoiss,CPoiss,BlocksNum,MaxIt)
    %Initialization
    AGamma = max(2^-52,rand(BlocksNum,BlocksNum));
    APoiss = max(2^-52,rand(BlocksNum,BlocksNum));

    for i = 1:MaxIt
        AGamma = max(2^-52, AGamma.*( ( CGamma*(((CGamma'*AGamma).^(-2)).*BGamma') )./( CGamma*((CGamma'*AGamma).^(-1)) ) ).^(0.5) );
        APoiss = max(2^-52, APoiss.*( ( CPoiss*(((CPoiss'*APoiss).^(-1)).*BPoiss') )./( repmat(sum(CPoiss,2),1,size(APoiss,2)) ) ) );
    end
    AGamma = AGamma';
    APoiss = APoiss'; 
end

function GeoMeans = ComputeGeoMeans(TestingSetOutput,InstancesNum,BlocksSize,BlocksNum)
    CSBlockSize = [0,cumsum(BlocksSize)];
    
    %Compute active indices 
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
end