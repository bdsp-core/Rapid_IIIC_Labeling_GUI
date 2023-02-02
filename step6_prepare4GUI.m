%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: Prepare input for GUI for patient with patient token "subject01"                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear

addpath('.\tools\preGUIcallbacks\')
eegDir   = '.\Data\processed\';       
qeegDir  = '.\Data\Spectrograms\';  
scoreDir = '.\Data\iiic\model_prediction\';

%% Define patient token - edit here to reflect your case
pat_token = 'subject01';

%% Define output folder 
trtDir = '.\Task\';
if exist(trtDir, 'dir')~=7
    mkdir(trtDir)
end

mainDir = [trtDir, pat_token, '\'];

disp('-------------------------------------------------------------------')
if exist(mainDir, 'dir')~=7
    disp(['step 0: make folders for ', pat_token])
    mkdir(mainDir)
    mkdir([mainDir, '\Data\EEG\']);
    mkdir([mainDir, '\Data\Labels\']);
    mkdir([mainDir, '\Data\Spectrograms\']);
    mkdir([mainDir, '\Data\Embedding\']);
end
dataDir_trt    = [mainDir, '\Data\EEG\'];
labelDir_trt   = [mainDir, '\Data\Labels\'];
specDir_trt    = [mainDir, '\Data\Spectrograms\'];
embedDir_trt   = [mainDir, '\Data\Embedding\'];

%% Get the list of files
files = struct2cell(dir([eegDir, '*.mat']))';
Files = strrep(files(:, 1), '.mat', '');

%% find EEG belong to that patient 
idx_selected = find(cellfun(@isempty, regexpi(Files, pat_token))==0);
Files = Files(idx_selected);

%% step 1: include EEG   
disp(['step 1: Add EEG for ', pat_token])

files = strrep(Files, '.mat', '');
fileMapping = cell(length(files), 2);
for k = 1:length(Files)
    file = files{k};     
    file1 = [pat_token, '_', repmat('0', 1, 2-length(num2str(k))), num2str(k)]; 
    fileMapping(k, :) = {file1, file};

    eegFile   = [file, '.mat'];   
    eegFile1  = [file1,'.mat'];  
    
    % copy file over
    disp(['  - copy file ', eegFile1])
    file_to_copy_src = [eegDir, eegFile];
    file_to_copy_trt = [dataDir_trt, eegFile1];
    copyfile(file_to_copy_src, file_to_copy_trt)
end  
disp('-------------------------------------------------------------------')

%% step 2: get labels and features
disp('step 2: Get SPaRCNet labels...')
labelDir_src = scoreDir;
for k = 1:size(files, 1)
    file = files{k};     
    file1 = [pat_token, '_', repmat('0', 1, 2-length(num2str(k))), num2str(k)]; 

    labelFile   = [file, '_score.mat'];   
    labelFile1  = [file1, '_score.mat'];   

    % copy file over
    disp(['  - copy file ', labelFile1]);
    file_to_copy_src = [labelDir_src, labelFile];
    file_to_copy_trt = [labelDir_trt, labelFile1];
    copyfile(file_to_copy_src, file_to_copy_trt);
end
disp('-------------------------------------------------------------------')

%% step 3: Get spectrograms 
disp('step 3: Get spectrograms...')
specDir_src = qeegDir;
for k =1:size(files, 1)
    file = files{k};     
    file1 = [pat_token, '_', repmat('0', 1, 2-length(num2str(k))), num2str(k)]; 

    specFile  = [file, '_spect.mat'];
    specFile1 = [file1, '_spec.mat'];

    % copy file over
    disp(['  - copy file ', specFile1])
    file_to_copy_src = [specDir_src, specFile];
    file_to_copy_trt = [specDir_trt, specFile1];   
    copyfile(file_to_copy_src, file_to_copy_trt)
end

%% step 3: Get PaCMAP
pacmapFile_out = [pat_token,'_pacmap.mat'];
pacmapFile_in  = [pat_token,'_feature.mat'];

disp('-------------------------------------------------------------------')
disp('step 4: Compute embedding map...');
disp([' - get normalized feature array ', pacmapFile_in]);

files = struct2cell(dir([dataDir_trt, '*.mat']))';
X_ = []; y_ = []; z_ = []; 
for k = 1:size(files, 1)
    file = strrep(files{k}, '.mat', '');
    scoreFile = [file, '_score.mat']; 

    disp(['  - add file ', file])
    tmp = load([labelDir_trt, scoreFile]);
    x = tmp.Y_model;

    [~, y] = max(x(:, [2:6, 1]), [], 2);
    y(y==6) = 0;

    X_ = [X_; x];                   % model score %
    y_ = [y_; y];                   % model prediction %
    z_ = [z_; k*ones(length(y), 1)];% EEG idx %
end
Y_model = X_;
y_model = y_;
idx_epoch = z_;

X = [X_, sum(X_(:, 2:4), 2)];
save([embedDir_trt,'\', pacmapFile_in], 'X', 'Y_model', 'y_model', 'idx_epoch')   

disp([' - compute pacmap ', pacmapFile_out]);
str = strrep(strrep(['"', embedDir_trt,'\', pacmapFile_in, ' "', embedDir_trt,'\',pacmapFile_out], '.mat', '"'),'\','/');   
commandStr = ['python ./Tools/computePaCMAP.py ', str];

% disp(commandStr)
system(commandStr);

% Export
tmp = load([embedDir_trt,'\', pacmapFile_out]);
Vxy = double(tmp.X_pacmap);
save([embedDir_trt,'\', pacmapFile_out], 'Vxy', 'Y_model', 'y_model', 'idx_epoch', 'fileMapping')   
disp('-------------------------------------------------------------------')

%%%% copy GUI %%%%
disp('step5: Copy GUI...')
gui_src = '.\Tools\CMGUI_Sequential_BoWspreading_v3.m';
gui_trt = [mainDir, 'CMGUI_Sequential_BoWspreading_v3.m'];
copyfile(gui_src, gui_trt)

disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('You are all set!')
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
 