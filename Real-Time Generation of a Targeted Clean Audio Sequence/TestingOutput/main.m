%DEEP NMF
%******************************************************************************************************
%% START
%MAKE SURE THE LATEST WEIGHTS HAVE BEEN SAVED BEFORE RUNNING AGAIN
clear;
close all;
clc;

fprintf('********************************************************************************\n');
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEEP NMF~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('********************************************************************************\n\n');

%*******************************************************************************************************
%% LOADING THE NET PARAMETERS AND INITIALIZING PROPAGATION OPTIONS
%Always save the net parameters using: save('Value of WorkingNetParametersName','NetParameters')
WorkingNetParametersName = 'NetParameters_Dataset3_PADNMF'; %Edit this to current working file name
 
fprintf('********************************************************************************\n');
fprintf('Loading net parameters contained in: %s\n', WorkingNetParametersName);
fprintf('********************************************************************************\n\n');

load(WorkingNetParametersName); %Creates variable 'NetParameters'

%Initializing propagation options (optional)
PropagationOptions = struct;
PropagationOptions.Epsilon = 10^-16;

%******************************************************************************************************
%% TESTING
%Always save the testing set using: save('Value of WorkingTestingSetName','TestingSet')
%Always save the testing set output using: save('Value of WorkingTestingSetOutputName','TestingSetOutput')
WorkingTestingSetName = 'TestingSetPiano_Dataset3'; %Edit this to current working file name
WorkingTestingSetOutputName = 'TestingSetSpindleOutput_Dataset3_PADNMF_'; %Edit this to current working file name
WorkingNetWeightsName = 'NetWeights_Dataset3_PADNMF'; %Edit this to current working file name

%-------------------------------------------------------------------------- 
fprintf('********************************************************************************\n');
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('********************************************************************************\n\n');

%Load the testing set 
%BEWARE: at the moment we are assuming the testing set to be a cell array
%of size (1,TestingInstancesNum), one cell dedicated to each testing instance. 
%Testing instances are mixture matrices.
fprintf('********************************************************************************\n');
fprintf('Loading testing set contained in: %s.\n', WorkingTestingSetName);
fprintf('********************************************************************************\n\n');

load(WorkingTestingSetName); %Creates variable 'TestingSet'

fprintf('********************************************************************************\n');
fprintf('Loading and using the net weights contained in: %s.\n', WorkingNetWeightsName);
fprintf('********************************************************************************\n\n');

load(WorkingNetWeightsName); %Creates variable 'NetWeights'

%Selecting the final dictionaries used for the reconstruction
RowNum = size(TestingSet{1},1);
WFinal = NetWeights{end}(end-RowNum+1:end,:);
WFinalSource = WFinal(end-RowNum+1:end,1:NetParameters.Ranks(1)); %Only the first source is reconstructed

%Detecting the number of instances
TestingInstancesNum = size(TestingSet,2);

%Creating a cell array of size (4,TestingInstancesNum) to contain the output,  
%one cell column dedicated to each output instance. In particular, the first row
%for each instance contains the reconstructed clean source (spectrogram),
%the second row for each instance contains the weights for an STFT weiner filter
%separating the clean source, the third row for each instance contains the
%reconstructed and 'wiener filtered' clean source (spectrogram) and the fourth row
%contains a cell with the complete history of H matrices.
TestingSetOutput = cell(4,TestingInstancesNum);

%Cycle over each testing instance
for TestingInstanceCounter = 1:TestingInstancesNum
    ShowProgress = fprintf('Currently processing testing instance: %d/%d.',TestingInstanceCounter,TestingInstancesNum);
    
    %Select instance mixture matrix
    MixMat = TestingSet{TestingInstanceCounter};
    %Propagate input
    HList = PropagateInput(MixMat,PropagationOptions,NetParameters,NetWeights);
    %Assemble output using the last coefficient matrix HList{end}
    TestingSetOutput{1,TestingInstanceCounter} = WFinalSource*HList{end}(1:NetParameters.Ranks(1),:);
    TestingSetOutput{2,TestingInstanceCounter} = TestingSetOutput{1,TestingInstanceCounter}./(WFinal*HList{end}+eps);
    TestingSetOutput{3,TestingInstanceCounter} = TestingSetOutput{2,TestingInstanceCounter}.*MixMat;
    for HCounter = 1:NetParameters.DiscriminativeLayers+1
        TestingSetOutput{4,TestingInstanceCounter}{HCounter} = HList{end-NetParameters.DiscriminativeLayers-1+HCounter};
    end
    fprintf(repmat('\b',1,ShowProgress));
end

%Save the output
save(WorkingTestingSetOutputName,'TestingSetOutput');

fprintf('********************************************************************************\n');
fprintf('Testing completed.\nThe output has been saved in %s and is ready to be analized.\n', WorkingTestingSetOutputName);
fprintf('********************************************************************************\n\n');

%Clear useless allocated variables
clear -regexp WorkingTestingSetName WorkingTestingSetOutputName WorkingNetWeightsName TestingSet NetWeights WFinal WFinalSource RowNum TestingInstancesNum TestingSetOutput TestingInstanceCounter ShowProgress MixMat HList HCounter;
%--------------------------------------------------------------------------

%******************************************************************************************************




