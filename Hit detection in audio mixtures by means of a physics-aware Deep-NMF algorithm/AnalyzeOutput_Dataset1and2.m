close all;
clear; 
clc;

%% ------------------------------------------------------CHOOSE OUTPUT TO REPRODUCE------------------------------------------------------------------ %%
%Configuration of the testing set
Config = 1; %1: Dataset1 - 2 macroblocks; 
            %2: Dataset2 - 5 macroblocks; 
            %3: Dataset2 - 5 macroblocks & 2 hits; 

if Config == 1
    %Dataset1 - 2 macroblocks
    TrainingMassesBlocks = {{'1',1},{'2',1},{'3',1},{'4',1},{'5',1},{'6',1},{'7',1},{'8',1},...
                    {'14',2},{'15',2},{'16',2},{'17',2},{'18',2},{'19',2},{'20',2},{'21',2},{'22',2},{'23',2}};  
    TrainingMasses = length(TrainingMassesBlocks);
    Ranks = 8*ones(1,TrainingMasses);
    CSRanks = [0,cumsum(Ranks)]; 
    Blocks = 2;
    
    TestingMassesBlocks = {{'1',1},{'2',1},{'3',1},{'4',1},{'5',1},{'6',1},{'7',1},{'8',1},...
                    {'14',2},{'15',2},{'16',2},{'17',2},{'18',2},{'19',2},{'20',2},{'21',2},{'22',2},{'23',2}};    
    TestingMasses = length(TestingMassesBlocks);           

    Strengths = {'Strong','Medium','Intermediate','Faint','Weak','VeryWeak','Absent'};
    TestingIntensities = length(Strengths);
    
    %Dataset directory
    DataDir = 'RawData/Dataset1DeepNMF/';

    %Output to show
    load('TestingOutput/TestingSetSpindleOutput_19ImpulsesPADNMF.mat'); %PADNMF
    %load('TestingOutput/TestingSetSpindleOutput_19ImpulsesDNMF.mat'); %DNMF

end

if Config == 2
    %Dataset2 - 5 macroblocks
    TrainingMassesBlocks = {{'2',1},{'3',1},...
                    {'8',2},{'9',2},...
                    {'14',3},{'15',3},...
                    {'20',4},{'21',4},...
                    {'26',5},{'28',5}}; 
    TrainingMasses = length(TrainingMassesBlocks);
    Ranks = 8*ones(1,TrainingMasses);
    CSRanks = [0,cumsum(Ranks)]; 
    Blocks = 5;
    
    TestingMassesBlocks = {{'2',1},{'3',1},...
                    {'8',2},{'9',2},...
                    {'14',3},{'15',3},...
                    {'20',4},{'21',4},...
                    {'26',5},{'28',5}};   
    TestingMasses = length(TestingMassesBlocks);           

    Strengths = {'Strong','Medium','Intermediate','Faint','Weak','VeryWeak','Absent'};
    TestingIntensities = length(Strengths);

    %Dataset directory
    DataDir = 'RawData/Dataset2DeepNMF/';

    %Output to show
    load('TestingOutput/TestingSetSpindleOutput_10ImpulsesPADNMF.mat'); %PADNMF
    
end

if Config == 3
    %Dataset2 - 5 macroblocks & 2 hits
    TrainingMassesBlocks = {{'2',1},{'3',1},...
                    {'8',2},{'9',2},...
                    {'14',3},{'15',3},...
                    {'20',4},{'21',4},...
                    {'26',5},{'28',5}}; 
    TrainingMasses = length(TrainingMassesBlocks);
    Ranks = 8*ones(1,TrainingMasses);
    CSRanks = [0,cumsum(Ranks)]; 
    Blocks = 5;
    
    TestingMassesBlocks = {{'2_8',1,2},...
                    {'8_14',2,3},...
                    {'14_20',3,4},...
                    {'20_26',4,5},...
                    {'2_26',5,1}};   
    TestingMasses = length(TestingMassesBlocks);

    Strengths = {'Medium','Absent'};
    TestingIntensities = length(Strengths);

    %Dataset directory
    DataDir = 'RawData/Dataset2DeepNMF/';

    %Output to show
    load('TestingOutput/TestingSetSpindleOutput_10ImpulsesPADNMFPairs.mat'); %PADNMF
    %load('TestingOutput/TestingSetSpindleOutput_10ImpulsesDNMFPairs.mat'); %DNMF
    
end


%% ------------------------------------------------------HIT DETECTION------------------------------------------------------------------ %%
BlockComponentsPosition = cell(1,Blocks);
for i = 1:TrainingMasses
    BlockComponentsPosition{TrainingMassesBlocks{i}{2}} = [BlockComponentsPosition{TrainingMassesBlocks{i}{2}},i];
end

%Output colors and symbols
colors = ['r','m','b','c','g']; %As many as (or more than) Blocks
symb = {'*','o','square','diamond','pentagram'}; %As many as (or more than) Blocks

%Load engine forcing for SNR calculation
Noise = importdata([DataDir,'EngineForcing.mat']);
Noise = Noise(0.4*44100:1.5*44100);

%For thresholding
NoImpulse = cell(1,Blocks);
Thresh = cell(1,Blocks);

%Will contain the geometric means for each testing instance
Coeffs = cell(2,Blocks);

StdNum = 1;
Off = 1e8;

ColLow = 0.1;
ColHigh = 0.7;

for i = 1:TestingMasses
    for j = 1:TestingIntensities
        ColList = linspace(ColLow,ColHigh,Blocks-1);
        
        if j < TestingIntensities
            load([DataDir, 'Impulse', TestingMassesBlocks{i}{1}, Strengths{j}, 'Forcing.mat']);
            SNR = snr(displ,Noise);
        end
        
        if Config ~= 3
            figure(i);
        else 
            figure(1);
        end
        
        for BlockInd = 1:Blocks
            BlockMassesInds = BlockComponentsPosition{BlockInd};
            for MassInd = 1:length(BlockMassesInds)
                Coeffs{1,BlockInd}(MassInd) = GeoMean(diag(TestingSetOutput{4,(i-1)*TestingIntensities+j}{end}(1+CSRanks(BlockMassesInds(MassInd)):end,1:Ranks(BlockMassesInds(MassInd)))));
            end
            
            Coeffs{2,BlockInd} = GeoMean(Coeffs{1,BlockInd});
            
            if j < TestingIntensities
                Mean = mean(10*log10(Off*Coeffs{1,BlockInd}));
                Std = std(10*log10(Off*Coeffs{1,BlockInd}));

                if (BlockInd == TestingMassesBlocks{i}{2})||(Config == 3 && BlockInd == TestingMassesBlocks{i}{3})
                   if BlockInd == 5
                       MS = 10;
                   else
                       MS = 9;
                   end
                   LW = 2; 
                   COL = [0,0,0];
                else
                   LW = 1;
                   MS = 7;
                   COL = ColList(1)*ones(1,3);
                   ColList = ColList(2:end);
                end
                
                plot(SNR,Mean,symb{BlockInd},'Color',COL,'MarkerSize',MS,'LineWidth',LW);
                hold on;
                if Config ~= 3
                    plot(SNR,Mean+StdNum*Std,'v','Color',COL,'MarkerSize',MS,'LineWidth',LW);
                    plot(SNR,Mean-StdNum*Std,'^','Color',COL,'MarkerSize',MS,'LineWidth',LW);
                end
                
            else
                NoImpulse{BlockInd}(i) = Coeffs{2,BlockInd};
            end
        end
    end
   
   xlabel('SNR (dB)','FontSize',18);
   ylabel('Geometric mean (dB)','FontSize',18);
end

%Computing the thresholds
for BlockInd = 1:Blocks
    Thresh{BlockInd} = GeoMean(NoImpulse{BlockInd});
end

%Plotting the thresholds
for i = 1:TestingMasses
   ColList = linspace(ColLow,ColHigh,Blocks-1);
   
   if Config ~= 3
        figure(i);
   else 
        figure(1);
   end

   for BlockInd = 1:Blocks
       
       if (BlockInd == TestingMassesBlocks{i}{2})
          LW = 2; 
          COL = [0,0,0];
       else
          LW = 1;
          COL = ColList(1)*ones(1,3);
          ColList = ColList(2:end);
       end
       yline(10*log10(Off*Thresh{BlockInd}),'-.','Color',COL,'LineWidth',LW);
   end
   grid on;
   grid minor;
end

%Geometric mean auxiliary function
function GM = GeoMean(v)
    mid = floor(length(v)/2);
    v1 = v(1:mid);
    v2 = v(mid+1:end);
    GM = ( prod(v1)^(1/length(v)) )*( prod(v2)^(1/length(v)) ); %to avoid underflow
end