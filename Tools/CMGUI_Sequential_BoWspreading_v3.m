function CMGUI_Sequential_BoWspreading_v3
    %-------------------------------------------------------------------------%
    % step 1: clustering into 50 clusters on model probability in K^2 distance
    % step 2: sequential inspection from the begining
    % 2-level spreading: BoW K^2 within 0.1 radius, then CP spreading
    % Accept the current entire window in spectrogram
    % All one class
    %-------------------------------------------------------------------------%
 
    close all
    clc
    warning('off','all')
   
    % Check the user crediential %
    user_ID = inputdlg('Please enter your initials','Hi');  
    rater = user_ID{1};
    targetDir = [pwd, '/Output_new/',rater,'/'];
    
    dataDir =  './Data/EEG/';
    embedDir = './Data/Embedding/';
    spectDir = './Data/Spectrograms/';  
    eegFiles = struct2cell(dir([dataDir, '*.mat']))';
    eegFiles = eegFiles(:, 1);

    f = figure('units','normalized','position',[0    0.0333    1.0000    0.9475]);
    set(f, 'MenuBar', 'none');
    set(f, 'ToolBar', 'none');
    set(f,'HitTest','off')
    set(f,'WindowButtonDownFcn',@clicks_Callback);
    set(f,'ButtonDownFcn',@clicks_Callback);
    set(f,'KeyPressFcn',@keys_Callback);
    
    f1 = figure('units','normalized','outerposition',[0 0 .5 1]);
    set(f1, 'MenuBar', 'none');
    set(f1, 'ToolBar', 'none');
    set(f1,'HitTest','off')
    set(f1,'KeyPressFcn',@keys_Callback);
    set(f1, 'Visible', 'off');
    
    % figure #2
    set(0,'CurrentFigure',f1)
    Ax_spec_bar1 = subplot('position',[.03  .05 .26 .02]);
    Ax_spec1   =  {subplot('position',[.03  .77 .26 .180]);    
                   subplot('position',[.03  .58 .26 .180]);  
                   subplot('position',[.03  .39 .26 .180]);   
                   subplot('position',[.03  .20 .26 .180]);   
                   subplot('position',[.03  .10 .26 .090]);   
                   subplot('position',[.03  .96 .26 .01])};   
    Ax_EEG1     = subplot('position',[.32  .1 .60  .85]); 
    Ax_EEG_bar1 = subplot('position',[.32 .05 .60 .02]);
    
    bg1 = uibuttongroup(f1,'Title','EEG patterns','units','normalized','Position', [.93 .7 .06 .25],'SelectionChangedFcn',@fcn_tmdecision);
    patterns = {'Seizure', 'LPD', 'GPD', 'LRDA', 'GRDA', 'Other'};
    shortcuts = {'1', '2', '3', '4', '5', '6'};
    nP = length(patterns);
    r1 = cell(nP, 1);
    ro1 = uicontrol(bg1,'Style','radiobutton','String','None','units','normalized','Position',[.1 1 .8 .1],'HandleVisibility','off');
    for ir = 1:nP
        r1{ir} = uicontrol(bg1,'Style','radiobutton','String',['[',shortcuts{ir},'] ', patterns{ir}, ],'units','normalized','Position',[.1 (-ir+nP+1)/(1+nP) .8 .08],'HandleVisibility','off');
    end
    goBack = uicontrol(f1,'style','pushbutton', 'units','normalized','position',[0.93    0.67    0.05   0.025], 'string','Go back','callback',@fcn_goBack);          

    % figure #1
    set(0,'CurrentFigure',f)
    Ax_spec_bar = subplot('position',[.03  .05 .26 .02]);
    Ax_spec   =  {subplot('position',[.03  .47 .26 .095]);  
                  subplot('position',[.03  .37 .26 .095]);  
                  subplot('position',[.03  .27 .26 .095]);   
                  subplot('position',[.03  .17 .26 .095]);   
                  subplot('position',[.03  .10 .26 .065]);   
                  subplot('position',[.03  .57 .26 .01])};   
    Ax_EEG     = subplot('position',[.32  .1 .6  .85]); 
    Ax_EEG_bar = subplot('position',[.32 .05 .6 .02]);
    Ax_visual     =  subplot('position',[.03 .66 .26 .29]); 
    Ax_visual_bar =  subplot('position',[.03 .63 .26 .02]);
    startx = uicontrol(f,'style','pushbutton','units','normalized',  'position',[0.4690    0.58    0.0520    0.05], 'string','Start','fontsize',15,'callback',@fcn_start);      
    done     = uicontrol(f,'style','pushbutton',         'units','normalized','position',[0.93    0.40-0.20    0.05   0.04], 'string','Done','fontsize',14,'callback',@fcn_done);     
    tm = uicontrol(f,'style','pushbutton',      'units','normalized','position',[0.93    0.67    0.05   0.025], 'string','Temp. match', 'fontsize',12, 'callback',@fcn_tm);          
    allSZ  = uicontrol(f,'style','pushbutton','units','normalized',  'position', [0.93    0.63    0.05   0.025], 'string','All SZ','fontsize',12,'callback',@fcn_allSZ);      
    allLPD  = uicontrol(f,'style','pushbutton','units','normalized',  'position', [0.93    0.60    0.05   0.025], 'string','All LPD','fontsize',12,'callback',@fcn_allLPD);    
    allGPD  = uicontrol(f,'style','pushbutton','units','normalized',  'position', [0.93    0.57    0.05   0.025], 'string','All GPD','fontsize',12,'callback',@fcn_allGPD);    
    allLRDA  = uicontrol(f,'style','pushbutton','units','normalized',  'position', [0.93    0.54    0.05   0.025], 'string','All LRDA','fontsize',12,'callback',@fcn_allLRDA);    
    allGRDA  = uicontrol(f,'style','pushbutton','units','normalized',  'position', [0.93    0.51    0.05   0.025], 'string','All GRDA','fontsize',12,'callback',@fcn_allGRDA);    
    allOther  = uicontrol(f,'style','pushbutton','units','normalized',  'position', [0.93    0.48    0.05   0.025], 'string','All Other','fontsize',12,'callback',@fcn_allOther);    
    bg = uibuttongroup(f,'Title','EEG patterns','units','normalized','Position', [.93 .7 .06 .25],'SelectionChangedFcn',@fcn_decision);
    nP = length(patterns);
    r = cell(nP, 1);
    ro = uicontrol(bg,'Style','radiobutton','String','None','units','normalized','Position',[.1 1 .8 .1],'HandleVisibility','off');
    for ir = 1:nP
        r{ir} = uicontrol(bg,'Style','radiobutton','String',['[',shortcuts{ir},'] ', patterns{ir}, ],'units','normalized','Position',[.1 (-ir+nP+1)/(1+nP) .8 .08],'HandleVisibility','off');
    end
    
    winSize  = 10;              
    winSizeH = uicontrol('Style', 'popup','String', {'60 min', '30 min', '10 min', '2 min'},'units','normalized','Position', [.26 .58 .03 .025], 'value', 3,'Callback', @fcn_winSize);    
    R_bow  = 0.1;              
    radiusSizeH = uicontrol('Style', 'popup', 'String', {'1.0', '0.8', '0.5', '0.3', '0.2', '0.1', '0.01'}, 'units','normalized', 'Position', [.3 .95 .03 .025], 'value', 6,  'Callback', @fcn_radiusSize);                  
    montage = 'L-Bipolar';
    Montage = uicontrol('Style', 'popup',  'String', {'L-Bipolar','Average','Monopolar'},  'units','normalized', 'Position', [.4 .95 0.045 0.0250], 'value', 1, 'Callback', @fcn_montageFlag);                     
    acceptNnext = uicontrol(f,'style','pushbutton',  'units','normalized','position',     [0.93    0.40    0.05   0.025], 'string', 'Accept all','fontsize',12, 'callback',@fcn_acceptNnext);   
    nextSPcenter = uicontrol(f,'style','pushbutton',      'units','normalized','position',[0.93    0.37    0.05   0.025], 'string', 'Next sample','fontsize',12,'callback',@fcn_goNextSPcenter);      
    progressStr = uicontrol(f,'style','text',      'units','normalized',  'position',     [0.93    0.30    0.07   0.025], 'string', 'Progress: -%','fontsize',12,'HorizontalAlignment','left');      

    set(Ax_EEG, 'Visible', 'off');    
    set(Ax_EEG_bar, 'Visible', 'off');    
    set(Ax_spec_bar, 'Visible', 'off');    
    for j = 1:length(Ax_spec)
        set(Ax_spec{j}, 'Visible', 'off'); 
    end 
    set(Ax_visual, 'Visible', 'off');
    set(Ax_visual_bar, 'Visible', 'off'); 
    set(bg, 'Visible', 'off'); 
    set(done, 'Visible', 'off'); 
    set(nextSPcenter, 'Visible', 'off'); 
    set(winSizeH, 'Visible', 'off');
    set(radiusSizeH, 'Visible', 'off');
    set(Montage, 'Visible', 'off');
    set(tm, 'Visible', 'off');
    set(acceptNnext, 'Visible', 'off');
    set(progressStr, 'Visible', 'off');
    set(allSZ, 'Visible', 'off'); set(allLPD, 'Visible', 'off'); set(allGPD, 'Visible', 'off'); set(allLRDA, 'Visible', 'off'); set(allGRDA, 'Visible', 'off');set(allOther, 'Visible', 'off');

    % parameters and variables
    colorMode = 'By pattern'; 
    colormap jet
    col = [-10 25];
    spatialRegs = {'LL', 'RL', 'LP', 'RP'};
    channel_withspace_bipolar  = {'Fp1-F7' 'F7-T3' 'T3-T5' 'T5-O1' '' 'Fp2-F8' 'F8-T4' 'T4-T6' 'T6-O2' '' 'Fp1-F3' 'F3-C3' 'C3-P3' 'P3-O1' '' 'Fp2-F4' 'F4-C4' 'C4-P4' 'P4-O2' '' 'Fz-Cz'  'Cz-Pz' '' 'EKG'};
    channel_withspace_average = {'Fp1-avg' 'F3-avg' 'C3-avg' 'P3-avg' 'F7-avg' 'T3-avg' 'T5-avg' 'O1-avg' '' 'Fz-avg' 'Cz-avg' 'Pz-avg' '' 'Fp2-avg' 'F4-avg' 'C4-avg' 'P4-avg' 'F8-avg' 'T4-avg' 'T6-avg' 'O2-avg' '' 'EKG'};
    channel_withspace_monopolar = {'Fp1' 'F3' 'C3' 'P3' 'F7' 'T3' 'T5' 'O1' '' 'Fz' 'Cz' 'Pz' '' 'Fp2' 'F4' 'C4' 'P4' 'F8' 'T4' 'T6' 'O2' '' 'EKG'}; 
    Fs = 200;
    kMediods = 50;
    K_bow = 50;
    thr_cp = .1;
    zScale = 1/150;   
    w = 14;            
    spec_step = 2;  
    tmFlag = 0;
    queryFlag = 1;
    doneFlag = [];
    Flags = [];  
    D_bow = [];
    BoW = [];
    bow_vec = [];
    y_tm = [];
    Y_model = [];
    DATA = [];
    Xvisual = [];    
    LUTtopN = [];
    LUTtopN_current = [];
    LUT_previous = []; 
    indTopN = [];
    idxReal = cell(0, 4);
    indSeen_   = [];   
    indUnseen_ = [];
    idxCPCs = [];
    idxCPs = [];
    isCPCs = [];
    isCPs = [];
    idx_mediods_cpc = [];
    idx_member2medoid_cpc = []; 
    ind_pairs = [];
    nSamples = [];
    jjj = [];
    jj = [];
    jjj1 = [];
    
    uiwait(f);
    while true
%         try
            set(done, 'Enable', 'on');
            set(nextSPcenter, 'Enable', 'on');
            set(winSizeH, 'Enable', 'on');
            set(radiusSizeH, 'Enable', 'on');
            set(Montage, 'Enable', 'on');
            set(acceptNnext, 'Enable', 'on');
            set(tm, 'Enable', 'on');
            set(allSZ, 'Enable', 'on');set(allLPD, 'Enable', 'on');set(allGPD, 'Enable', 'on');set(allLRDA, 'Enable', 'on');set(allGRDA, 'Enable', 'on');set(allOther, 'Enable', 'on');         
            for j = 1:nP
                set(r{j}, 'Enable', 'on')
            end
            set(ro, 'Enable', 'on');
            set(0,'CurrentFigure',f)
            pp = sum(doneFlag==1)/length(doneFlag);            
            set(progressStr, 'String', ['Progress: ', num2str(round(1E4*pp)/100), '%']);

            iEpoch = LUT_previous{jjj, 1};
            ii     = LUT_previous{jjj, 2};
            if isempty(find(indTopN==jjj, 1))
                y_hat = [];
            else
                y_hat = LUTtopN_current{find(indTopN==jjj, 1), 6};
            end
   
            data   = DATA{iEpoch, 2};
            Sdata  = DATA{iEpoch, 3}{1}(:, 2);
            stimes = DATA{iEpoch, 3}{2}; 
            sfreqs = DATA{iEpoch, 3}{3};
            totalPower = DATA{iEpoch, 4}{2}; 
            cpds = DATA{iEpoch, 4}{1};
 
            idx_local = find(cell2mat(LUT_previous(:, 1)) == iEpoch);    
            if queryFlag == 1  
                [~, Y_clustercpdSmooth] = fcn_colorClusterCPDsmooth;  
                LUT_previous(:, 6)= Y_clustercpdSmooth;
                labels = Y_clustercpdSmooth(idx_local);
            else
                Y_bowcpdSmooth = fcn_labelBoWCPDsmooth;
                LUT_previous(:, 6)= Y_bowcpdSmooth;
                labels = Y_bowcpdSmooth(idx_local);
                disp('Sequential model')
            end

            idx_current_local = find((doneFlag==1) & (cell2mat(LUTtopN_current(:, 1)) == iEpoch));  
            ii_done =  cell2mat(LUTtopN_current(idx_current_local, 2)); 
            labels(ii_done) = LUTtopN_current(idx_current_local, 6);

            isCPs_local  = cell2mat(LUT_previous(idx_local, 3));              
            isCPCs_local = cell2mat(LUT_previous(idx_local, 4));  

            idxCPCs_local = find(isCPCs_local==1);  
            idxCPs_local  = find(isCPs_local~=0);  
            for j = 1:length(idxCPCs_local)
                idx1 = idxCPs_local(find(idxCPs_local<=idxCPCs_local(j), 1, 'last')); 
                idx2 = -1+idxCPs_local(find(idxCPs_local>idxCPCs_local(j), 1, 'first')); 
                scr_cpc = labels{idxCPCs_local(j)};
                labels(idx1:idx2) = repmat({scr_cpc}, (idx2-idx1+1), 1);
            end
            labels(end) = labels(end-1);

            barColors = repmat([.8 .8 .8], length(labels), 1) ;
            Colors6 = flipud(jet(7));Colors6 = Colors6(1:6, :);
            ystr = categorical(labels); 
            for j = 1:length(patterns)
                idx_j = find(ystr == patterns{j});
                barColors(idx_j, :) = repmat(Colors6(j,:), length(idx_j), 1);
            end
            tc  = 2*ii-1;         
            itc = tc*Fs+1;       
            iStart = max(1, itc-Fs*w/2+1);          
            iEnd   = min(size(data, 2), itc+Fs*w/2);          
            seg = fcn_parseData(data(:, iStart:iEnd), itc, iStart, iEnd,  size(data, 1), w*Fs, size(data, 2));
             
            barColors_eeg_ = barColors(max(ii-(w/2-1)/2, 1) : min(length(labels), ii+(w/2-1)/2), :);
            barColors_eeg = ones(w/2, 3);
            if ii-(w/2-1)/2<1
                barColors_eeg(w/2-size(barColors_eeg_, 1)+1 : w/2, :) = barColors_eeg_;
            elseif ii+(w/2-1)/2> length(labels)
                barColors_eeg(1: size(barColors_eeg_, 1), :) = barColors_eeg_;
            else
                barColors_eeg = barColors_eeg_;
            end

            W = winSize*60;
            iStart_f = max(1, ii-W/(2*spec_step));  
            iEnd_f   = min(length(stimes), ii+W/(2*spec_step)-1);
            Sdata_ = cell(4, 1);
            for j = 1:4
                S = Sdata{j};      
                seg_f = S(:, iStart_f:iEnd_f);
                if j == 4
                    barColors_spec_ = barColors(iStart_f:iEnd_f, :);
                    total_power_ = totalPower(iStart_f:iEnd_f);
                end
                if size(seg_f, 2) == (W/spec_step)
                    spec = seg_f;
                    if j == 4
                        barColors_spec = barColors_spec_;
                        total_power = total_power_;
                    end
                else  
                    d = (W/spec_step) - size(seg_f, 2);
                    if iStart_f == 1 && iEnd_f ~= length(stimes)
                        spec = [eps+zeros(size(seg_f, 1), d), seg_f];
                        if j == 4
                            barColors_spec = [ones(d, 3); barColors_spec_];
                            total_power = [ones(1, d), total_power_];
                        end
                    elseif iEnd_f == length(stimes) && iStart_f ~= 1 
                        spec = [seg_f, eps+zeros(size(seg_f, 1), d)];
                        if j == 4
                            barColors_spec = [barColors_spec_; ones(d, 3)];
                            total_power = [total_power_, ones(1, d)];
                        end
                    elseif iEnd_f == length(stimes) && iStart_f == 1 
                        a = 1-(ii-W/(2*spec_step));
                        b = ii+ W/(2*spec_step) - length(stimes)-1;
                        spec = [eps+zeros(size(seg_f, 1), a), seg_f, eps+zeros(size(seg_f, 1), b)];
                        if j == 4
                            barColors_spec = [ones(a, 3); barColors_spec_; ones(b, 3)];
                            total_power = [ones(1, a), total_power_, ones(1, b)];
                        end
                    end
                end
                Sdata_{j, 1} = spec;
            end
            Sdata = Sdata_;

            timeStampo = DATA{iEpoch, 1}{2};            
            eeg = seg(1:19,:);
            ekg = seg(20, :);          
            tto = tc*Fs - (w/2)*Fs+1;   
            tt1 = tc*Fs + (w/2)*Fs;    
            tt = tto:tt1;
            gap = NaN(1, size(eeg, 2));
            switch montage
                case 'L-Bipolar'
                    seg = fcn_bipolar(eeg);
                    seg_disp = [seg(1:4,:); gap; seg(5:8,:); gap; seg(9:12,:); gap; seg(13:16,:); gap; seg(17:18,:); gap; ekg];
                    channel_withspace = channel_withspace_bipolar;
                case 'Average' 
                    seg = eeg - repmat(mean(eeg, 1), size(eeg, 1), 1);
                    seg_disp = [seg(1:8,:); gap; seg(9:11,:); gap; seg(12:19,:); gap; ekg];
                    channel_withspace = channel_withspace_average;
                case 'Monopolar'
                    seg =  eeg;
                    seg_disp = [seg(1:8,:); gap; seg(9:11,:); gap; seg(12:19,:); gap; ekg];
                    channel_withspace = channel_withspace_monopolar;
            end
            M = size(seg_disp, 1);
            DCoff = repmat(flipud((1:M)'), 1, size(seg_disp, 2));           
            timeStamps = datestr(timeStampo + seconds(round(tto/Fs):2:round(tt1/Fs)), 'hh:MM:ss');

            set(f,'CurrentAxes',Ax_EEG_bar); cla(Ax_EEG_bar)
            colors_ = reshape(barColors_eeg , 1, size(barColors_eeg, 1), 3); 
            image(Ax_EEG_bar,round(tt(1)/Fs)+(1:2:14), 1, colors_);
            set(Ax_EEG_bar,'ytick',[],'xtick',[])
            box(Ax_EEG_bar, 'on')

            set(f,'CurrentAxes',Ax_EEG); cla(Ax_EEG)
            hold(Ax_EEG, 'on')
                for j = 1:w
                    tj = tto + Fs*(j)-1;
                    line([tj  tj],  [0 M+1], 'linestyle', '--', 'color', [.5 .5 .5])
                end
                plot(Ax_EEG, tt, zScale*seg_disp(1:end-1,:)+DCoff(1:end-1,:),'k');
                ekg_ = seg_disp(end,:);
                ekg_ = (ekg_-mean(ekg_))/(eps+std(ekg_));

                plot(Ax_EEG, tt, .2*ekg_+DCoff(end,:),'m');
                box(Ax_EEG, 'on')
                set(Ax_EEG, 'ytick',1:M,'yticklabel',fliplr(channel_withspace),'box','on', 'ylim', [0 M+1], 'xlim', [tto-1 tt1], 'xtick',(tto-1:2*Fs:tt1),'xticklabel',timeStamps)

                dt = tt1-tto+1; a = round(dt*7.5/10);
                xa1 = tto+[a a+Fs-1]; ya1 = [2 2];
                xa2 = tto+[a a]; ya2 = ya1+[0 100*zScale];

                text(Ax_EEG, xa1(1)-Fs/20,  mean(ya2), '100\muV','Color', 'b','FontSize',15, 'horizontalalignment', 'right', 'verticalalignment', 'middle');
                text(Ax_EEG, mean(xa1), 1.95, '1 sec','Color', 'b','FontSize',15, 'horizontalalignment', 'center', 'verticalalignment', 'top');        
                line(Ax_EEG, xa1, ya1, 'LineWidth', 1, 'Color','b');
                line(Ax_EEG, xa2, ya2, 'LineWidth', 1, 'Color','b');
                if ~isempty(y_hat)
                    set(bg,'SelectedObject',r{ismember(patterns, y_hat)})
                else  
                    set(bg,'SelectedObject',ro)
                end
            hold(Ax_EEG, 'off')

            to = (tc)-(winSize*60/2);  
            t1 = (tc)+(winSize*60/2)-1;
            stepSize = 2;
            S_x = to:stepSize:t1;   
            S_y = sfreqs;
            for j = 1:4
                set(f,'CurrentAxes',Ax_spec{j}); cla(Ax_spec{j})
                spec = Sdata{j};              
                hold(Ax_spec{j}, 'on')
                    imagesc(Ax_spec{j}, S_x, S_y, pow2db(spec), col);  axis(Ax_spec{j},'xy'); 
                    box(Ax_spec{j}, 'on')
                    set(Ax_spec{j},'xtick',[], 'ylim', [min(S_y), max(S_y)])   
                    plot(Ax_spec{j}, [tc  tc], [S_y(1) S_y(end)], 'k--','linewidth',1);
                hold(Ax_spec{j}, 'off')
                xx = get(Ax_spec{j}, 'yticklabel'); xx{end} = spatialRegs{j};
                xlim([to, t1+1])
                set(Ax_spec{j}, 'yticklabel',xx)
                ylabel(Ax_spec{j}, 'Freq (Hz)')
            end

            set(f,'CurrentAxes',Ax_spec{5}); cla(Ax_spec{5})
            hold(Ax_spec{5}, 'on')  
                if length(S_x)~=length(total_power)
                    keyboard
                end
                plot(Ax_spec{5}, S_x, total_power, 'g-', 'linewidth', 1)
                cpds_tmp_ = stimes(cpds)-1;
                idx_tmp = find(cpds_tmp_>=to & cpds_tmp_<=t1);
                cpds_tmp = cpds_tmp_(max(1, min(idx_tmp)-1) : min(length(cpds_tmp_), max(idx_tmp)+1));
                for j = 1:length(cpds_tmp)-1
                    aa_ = cpds_tmp(j);
                    bb_ = cpds_tmp(j+1);
                    plot(Ax_spec{5}, [aa_, aa_], [-10 25], 'm-', 'linewidth', .5)
                    cpd_mean = mean(totalPower((aa_/2+1):(bb_/2+1)));
                    plot(Ax_spec{5}, [aa_, bb_], [cpd_mean, cpd_mean], 'b-', 'linewidth', .5)
                end
                plot(Ax_spec{5}, [tc  tc], [-10 25], 'k--','linewidth',1);
                tt = timeStampo + seconds(to:(winSize/5)*60:t1+1);
                set(Ax_spec{5}, 'xtick', to:(winSize/5)*60:t1+1, 'xticklabel', datestr(tt, 'hh:MM:ss'),...
                    'xlim',[to, t1+1],...
                    'ylim', [-10 25],...
                    'ytick', [], 'yticklabel', []);
                box(Ax_spec{5}, 'on')    
            hold(Ax_spec{5}, 'off')    

            set(f,'CurrentAxes',Ax_spec{end});
            cla(Ax_spec{end})
            plot(Ax_spec{end}, tc, 0,'rv','markersize', 8,'MarkerFaceColor','r'); 
            axis(Ax_spec{end}, 'off');
            set(Ax_spec{end},'xtick',[],'xlim',get(Ax_spec{end-1},'xlim'), 'ylim', [0 .5]) 

            set(f,'CurrentAxes',Ax_spec_bar);cla(Ax_spec_bar)
            colors = reshape(barColors_spec , 1, size(barColors_spec, 1), 3); 
            image(Ax_spec_bar, S_x, 1, colors);
            set(Ax_spec_bar,'xtick',[],'ytick',[])
            xlim([to, t1+1])
            box(Ax_spec_bar, 'on')

            set(f,'CurrentAxes', Ax_visual_bar); 
            colors = reshape(Colors6 , 1, size(Colors6 , 1), 3); 
            imagesc(Ax_visual_bar, 1:6,1,colors)
            set(Ax_visual_bar, 'xtick', 1:6, 'xticklabel', {'Seizure' 'LPD' 'GPD' 'LRDA' 'GRDA' 'Other'}, 'ytick', [])

            set(f,'CurrentAxes', Ax_visual); 
            Vxy = Xvisual;
            lut = LUT_previous; lut(indTopN(doneFlag==1), 6) = LUTtopN_current(doneFlag==1, 6);       
            if queryFlag == 1
                colors = fcn_colorClusterCPDsmooth;  
            else
                colors = fcn_pattern2color(lut, colorMode);   
            end           
            hold(Ax_visual, 'on'); cla(Ax_visual);
                ss = scatter(Ax_visual, Vxy(:, 1), Vxy(:, 2), 30, colors, 'filled'); 
                alpha(ss,.5)
                scatter(Ax_visual, Vxy(indSeen_, 1),  Vxy(indSeen_, 2), 30, colors(indSeen_, :), 'filled'); 
                plot(Ax_visual, Vxy(jjj, 1), Vxy(jjj, 2), 'mx','markersize', 20, 'linewidth', 2);
                plot(Ax_visual, Vxy(jjj, 1), Vxy(jjj, 2), 'mo','markersize', 20, 'linewidth', 2);
                scatter(Ax_visual, Vxy(jjj, 1), Vxy(jjj, 2), 30, colors(jjj,:), 'filled');
                idx_done = indTopN(doneFlag==1);  
                scatter(Ax_visual, Vxy(idx_done, 1), Vxy(idx_done, 2), 30, colors(idx_done, :), 'v','filled', 'MarkerEdgeColor','k');     
                axis(Ax_visual, 'off'); axis(Ax_visual, 'equal');
                set(Ax_visual, 'xlim', [min(Vxy(:, 1)), max(Vxy(:, 1))], 'ylim', [min(Vxy(:, 2)), max(Vxy(:, 2))])        
                if ~isempty(find(indTopN == jjj, 1))
                    idx_c = jjj;
                else
                    idx1_ = idxCPs(find(idxCPs<=jjj, 1, 'last')); 
                    idx2_ = idxCPs(find(idxCPs>jjj, 1, 'first')); 
                    if isempty(idx2_)  
                        idx_c = jjj;
                    else
                        idx_c = floor((idx1_+idx2_)/2); 
                    end
                end    
                if queryFlag == 1 
                    if ~isempty(find(indTopN == idx_c, 1))
                        idx_m = idxCPCs(idx_member2medoid_cpc== idx_mediods_cpc(indTopN == idx_c)); 
                        title(['Mode 1: Cluster #',  num2str(find((indTopN == idx_c))), '/50 with ', num2str(length(idx_m)), ' memebers'])    
                    else
                        title('Mode 1: non-center')    
                    end
                else
                    dd = D_bow(idxCPCs==idx_c, :);
                    idx_m = setdiff(idxCPCs(dd<R_bow), idx_done); 
                    title(['Model 2: SP center #',  num2str(find((indTopN == idx_c))), ' with ', num2str(length(find(dd<R_bow))), ' memebers'])
                end
                scatter(Ax_visual, Vxy(idx_m, 1), Vxy(idx_m, 2), 30, colors(idx_c, :), 'o','filled', 'MarkerEdgeColor','k'); 
                scatter(Ax_visual, Vxy(idx_c, 1), Vxy(idx_c, 2), 50, colors(idx_c, :), 's','filled', 'MarkerEdgeColor','k'); 
            hold(Ax_visual, 'off')
            
            uiwait(f);
            
%         catch err
%             fcn_done; 
%         end
    end
  
    function fcn_start(varargin) 
        set(startx,'Visible','off','Enable', 'off');
        loading = uicontrol(f,'style','text','units','normalized','position',[0.4690    0.5000+.08    0.0520    0.050],'string','Loading...','FontSize',15);
        drawnow

        DATA = cell(length(eegFiles), 4);
        Flags = [];
        
        file_tmp = struct2cell(dir([embedDir, '*pacmap.mat']))';
        tmp = load([embedDir, file_tmp{1}]);
        idx_epoch = tmp.idx_epoch;
        Xvisual = tmp.Vxy;
        Y_model = tmp.Y_model;

        for iFile = 1:length(eegFiles)
            dataFile = eegFiles{iFile};
            spectFile = strrep(dataFile, '.mat', '_spec.mat');
            tmp = load([dataDir, dataFile]);
            data = tmp.data;  
            data(isnan(data)) = 9999;
            startTime = tmp.startTime;          
            if ~isdatetime(startTime)
                startTime = datetime(startTime);
            end
            DATA{iFile, 1} = {strrep(dataFile, '.mat', ''), startTime};
            fc = 60;
            [B1, A1] = butter(3, [.5, 40]/(Fs/2));
            [B2, A2] = butter(3, [fc-2.5, fc+2.5]/(Fs/2), 'stop');            
            data = filtfilt(B1, A1, data')';
            data = filtfilt(B2, A2, data')';
            thr = 300;
            data(data<=-thr) = -thr;
            data(data>= thr) =  thr;
            DATA{iFile, 2} = data;
            tmp = load([spectDir, spectFile]);
            sdata_tmp = tmp.Sdata;
            for ireg = 1:4
                s = sdata_tmp{ireg, 2};
                sdata_tmp{ireg, 2} = [eps+zeros(size(s, 1), 1), s]; 
            end
            nn = size(sdata_tmp{ireg, 2}, 2);
            mm = sum(idx_epoch==iFile);
            if nn<=mm
                dd = mm-nn;
                for ireg = 1:4
                    s = sdata_tmp{ireg, 2};
                    sdata_tmp{ireg, 2} = [s, eps+zeros(size(s, 1), dd)];  
                end
            else
                for ireg = 1:4
                    s = sdata_tmp{ireg, 2};
                    sdata_tmp{ireg, 2} = s(:, 1:mm);  
                end
            end
            stime_tmp = 1+(0:2:(mm-1)*2);
            DATA{iFile, 3} = {sdata_tmp, stime_tmp, tmp.sfreqs};
            [cpd, p, cps, cpcs] = fcn_cpd(DATA{iFile, 3}{1}, mm);
            DATA{iFile, 4} = {cpd, p};
            Flags = [Flags; [cps cpcs]];
        end

        if exist([targetDir, 'RaterScores_final.mat'], 'file')==2 
            tmp_previous = load([targetDir, 'RaterScores_final.mat']);
            LUT_previous = tmp_previous.LUT;   
            LUT_previous = [LUT_previous, LUT_previous(:, end)];           
            indSeen_ = tmp_previous.idxSeen;           
            Y_model = tmp_previous.Y_model;
            idxReal = tmp_previous.idxReal; 
            thr_cp = tmp_previous.thr_cp;            
            for i = 1:size(DATA, 1)
                y_cpd = cell2mat(LUT_previous(cell2mat(LUT_previous(:, 1))==i, 3));
                DATA{iFile, 4}{1} = find(y_cpd~=0)';
            end
        else   
            if exist([targetDir, 'RaterScores_initial.mat'], 'file') ~= 2                
                file_tmp = struct2cell(dir([embedDir, '*pacmap.mat']))';
                tmp = load([embedDir, file_tmp{1}], 'y_model', 'idx_epoch');
                y_ =  tmp.y_model;  
                y_(y_==0) = 6;
                [~, Ystr] = fcn_label2number(y_, []);
                z_ = tmp.idx_epoch;  
                iz = []; z = unique(z_);
                for i_ = 1:length(z)
                    n = length(find(z_ ==z(i_)));
                    iz = [iz, 1:n];
                end             
                LUT = [num2cell([z_, iz',  Flags]), Ystr, Ystr];   
                LUT_previous = LUT;           
                if exist(targetDir, 'dir') ~= 7  
                    mkdir(targetDir)
                end
                save([targetDir,'RaterScores_initial.mat'], 'LUT', 'Y_model')               
            else
                tmp = load([targetDir,'RaterScores_initial.mat'], 'LUT', 'Y_model');
                LUT_previous = tmp.LUT; 
            end     
        end
        nSamples = size(LUT_previous, 1);
        isCPCs   = cell2mat(LUT_previous(:, 4)); 
        idxCPCs  = find(isCPCs==1);
        isCPs    = cell2mat(LUT_previous(:, 3)); 
        idxCPs   = find(isCPs~=0);       
        indUnseen_ = setdiff(idxCPCs, indSeen_);
        if isempty(indUnseen_)
            choice = questdlg('Wow! Seems you have labeled all the segments!','Warning','Ok, I am done.','Ok, I am done.');
            switch choice
                case 'Ok, I am done.'
                    fcn_done();
            end
            close all
            delete(f)
        else
            ind_pairs = cell2mat(LUT_previous(:, [1 2]));
            BoW = fcn_getBoW(DATA, idx_epoch, K_bow);            
            bow_vec = NaN(length(idxCPCs), K_bow);
            for m = 1:length(idxCPCs)
                idx1_ = idxCPs(find(idxCPs<=idxCPCs(m), 1, 'last'));  
                idx2_ = -1+idxCPs(find(idxCPs>idxCPCs(m), 1, 'first'));      
                bow_vec(m, :) = hist(BoW(idx1_:idx2_), 1:K_bow);
            end
            bow_vec = bow_vec./sum(bow_vec, 2);
            D_bow = fcn_distChiSq(bow_vec, bow_vec);  

            % Query mode 1 - step 1 clusters %
            queryFlag = 1;
            X = Y_model(idxCPCs, :);
            kMediods = min(length(idxCPCs), kMediods);
            [idx_mediods_cpc, idx_member2medoid_cpc] =  fcn_kmediodDistChiSq(X, kMediods);             
            indTopN = idxCPCs(idx_mediods_cpc);  

            % Get ready for re-lable %
            LUTtopN = LUT_previous(indTopN, :);
            N = length(indTopN);                                
            LUTtopN_current = cell(N, 6);                       
            LUTtopN_current(:, 1:5) = LUTtopN(:, [1:4, 6]);    
            doneFlag = NaN(N, 1);
            doneFlag(ismember(indTopN, indSeen_)) = 1;
            jj = find(isnan(doneFlag), 1);
            jjj = indTopN(jj);
            
            delete(loading)
            set(Ax_EEG, 'Visible', 'on');  
            set(Ax_EEG_bar, 'Visible', 'on');    
            set(Ax_spec_bar, 'Visible', 'on'); 
            for iAx = 1:length(Ax_spec)
                set(Ax_spec{iAx}, 'Visible', 'on'); 
            end
            set(Ax_visual, 'Visible', 'on'); 
            set(Ax_visual_bar, 'Visible', 'on'); 
            set(bg,     'Visible', 'on'); 
            set(winSizeH, 'Visible', 'on');
            set(radiusSizeH, 'Visible', 'on');
            set(Montage, 'Visible', 'on');
            set(progressStr, 'Visible', 'on');
            if queryFlag==1
                set(done, 'Visible', 'off');
                set(nextSPcenter, 'Visible', 'off');
                set(acceptNnext, 'Visible', 'off');
                set(tm, 'Visible', 'off');
                set(allSZ, 'Visible', 'off');
                set(allLPD, 'Visible', 'off');
                set(allGPD, 'Visible', 'off');
                set(allLRDA, 'Visible', 'off');
                set(allGRDA, 'Visible', 'off');
                set(allOther, 'Visible', 'off');
            else
                set(done, 'Visible', 'on');
                set(nextSPcenter, 'Visible', 'on');
                set(acceptNnext, 'Visible', 'on');
                set(tm, 'Visible', 'on');
                set(allSZ, 'Visible', 'on');
                set(allLPD, 'Visible', 'on');
                set(allGPD, 'Visible', 'on');
                set(allLRDA, 'Visible', 'on');
                set(allGRDA, 'Visible', 'on');
                set(allOther, 'Visible', 'on');
            end
            uiresume(f);
        end
    end

    function dataBipolar = fcn_bipolar(data)
        %labels = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};
        dataBipolar( 1,:) = data( 1,:) - data( 5,:);  % Fp1-F7
        dataBipolar( 2,:) = data( 5,:) - data( 6,:);  % F7-T3
        dataBipolar( 3,:) = data( 6,:) - data( 7,:);  % T3-T5
        dataBipolar( 4,:) = data( 7,:) - data( 8,:);  % T5-O1

        dataBipolar( 5,:) = data(12,:) - data(16,:);  % Fp2-F8
        dataBipolar( 6,:) = data(16,:) - data(17,:);  % F8-T4
        dataBipolar( 7,:) = data(17,:) - data(18,:);  % T4-T6
        dataBipolar( 8,:) = data(18,:) - data(19,:);  % T6-O2
        
        dataBipolar( 9,:) = data( 1,:) - data( 2,:);  % Fp1-F3
        dataBipolar(10,:) = data( 2,:) - data( 3,:);  % F3-C3
        dataBipolar(11,:) = data( 3,:) - data( 4,:);  % C3-P3
        dataBipolar(12,:) = data( 4,:) - data( 8,:);  % P3-O1

        dataBipolar(13,:) = data(12,:) - data(13,:);  % Fp2-F4
        dataBipolar(14,:) = data(13,:) - data(14,:);  % F4-C4
        dataBipolar(15,:) = data(14,:) - data(15,:);  % C4-P4
        dataBipolar(16,:) = data(15,:) - data(19,:);  % P4-O2

        dataBipolar(17,:) = data( 9,:) - data(10,:);  % Fz-Cz
        dataBipolar(18,:) = data(10,:) - data(11,:);  % Cz-Pz
    end

    function keys_Callback(f,varargin)
        k = get(f,'CurrentKey');
        disp(k)
        switch k
            case 'rightarrow'
                if tmFlag==1  
                    tmFlag = 0;  set(f1, 'Visible', 'off');          
                else
                    if queryFlag == 1
                        jj = find(indTopN==jjj);
                        jj = min(length(indTopN), jj+1);
                        jjj = indTopN(jj);                      
                    else
                        jjj = min(nSamples, jjj+1);                   
                    end                    
                end
            case 'leftarrow'
                if tmFlag==1 
                    tmFlag = 0;  set(f1, 'Visible', 'off');       
                else
                    if queryFlag == 1
                        jj = find(indTopN==jjj);
                        jj = max(1, jj-1);
                        jjj = indTopN(jj);                  
                    else
                        jjj = max(1, jjj-1);                
                    end
                end
            case 'uparrow'
                zScale =  zScale*1.5;             
            case 'downarrow'
                zScale =  zScale/1.5;            
            case {'1', 'numpad1'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{1})
                    fcn_tmdecisionkeys(1);
                else
                    set(bg,'SelectedObject',r{1})
                    fcn_decisionkeys(1);
                end
            case {'2', 'numpad2'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{2})
                    fcn_tmdecisionkeys(1);
                else
                    set(bg,'SelectedObject',r{2})
                    fcn_decisionkeys(2);
                end
            case {'3', 'numpad3'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{3})
                    fcn_tmdecisionkeys(1);
                else
                    set(bg,'SelectedObject',r{3})
                    fcn_decisionkeys(3);
                end      
            case {'4', 'numpad4'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{4})
                    fcn_tmdecisionkeys(1);
                else
                    set(bg,'SelectedObject',r{4})
                    fcn_decisionkeys(4);
                end
            case {'5', 'numpad5'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{5})
                    fcn_tmdecisionkeys(1);
                else
                    set(bg,'SelectedObject',r{5})
                    fcn_decisionkeys(5);
                end  
            case {'6', 'numpad6'}
                if tmFlag==1
                    set(bg1,'SelectedObject',r1{6})
                    fcn_tmdecisionkeys(1);
                else
                    set(bg,'SelectedObject',r{6})
                    fcn_decisionkeys(6);
                end     
            case {'space'}
                if tmFlag==1
                    tmFlag = 0; set(f1, 'Visible', 'off');
                else
                    if queryFlag == 1
                        fcn_nextCluster
                    else
                        fcn_nextSPcenter;
                    end
                end
            otherwise
                disp(['Dude! Invalid Key ', k])
        end
        uiresume(f);
    end

    function fcn_decisionkeys(ind_)
        p = patterns{ind_};   
        idxReal = [idxReal;{jjj, p, datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}];       
        idx1_ = idxCPs(find(idxCPs<=jjj, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>jjj, 1, 'first'));  
        if isempty(idx2_)  
            tmp_c = jjj;
        else
            tmp_c = floor((idx1_+idx2_)/2); 
        end
        if isempty(find(indTopN == tmp_c, 1))         
            keyboard
        else
            jj_c = find(indTopN == tmp_c); 
        end
        LUTtopN_current{jj_c, 6} = p;
        doneFlag(jj_c) = 1;    
        if tmFlag==-1
            fcn_tm();
        else
            if queryFlag == 1
                fcn_nextCluster
            else
                fcn_nextSPcenter;
            end
        end 
    end

    function fcn_decision(source, event)
        for i = 1:nP
            set(r{i}, 'Enable', 'off');
        end
        set(ro, 'Enable', 'off');
        drawnow;       
        p = event.NewValue.String;
        ip = regexp(p,']');
        p = p(ip+2:end);      
        idxReal = [idxReal; {jjj, p, datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}];
        idx1_ = idxCPs(find(idxCPs<=jjj, 1, 'last'));  
        idx2_ = idxCPs(find(idxCPs>jjj, 1, 'first'));  
        if isempty(idx2_)  
            tmp_c = jjj;
        else
            tmp_c = floor((idx1_+idx2_)/2); 
        end
        if isempty(find(indTopN == tmp_c, 1))  
            keyboard 
        else
            jj_c = find(indTopN == tmp_c);  
        end

        LUTtopN_current{jj_c, 6} = p;
        doneFlag(jj_c) = 1;
        
        if tmFlag==-1
            fcn_tm();
        else
            if queryFlag == 1
                fcn_nextCluster
            else
                fcn_nextSPcenter;
            end
        end      
    end

    function c = fcn_pattern2color(LUT, colorMode)
        c = .8*ones(size(LUT, 1), 3);
        switch colorMode
            case 'By pattern'
                Cs = flipud(jet(7));
                
                scr = categorical(LUT(:, 6));  
                for k = 1:length(patterns)
                    ind = find(scr == patterns{k});
                    c(ind, :) = repmat(Cs(k,:), length(ind), 1);
                end
        end
    end

    function clicks_Callback(varargin)
        set(bg,'SelectedObject',ro)
        click_type = get(f,'selectiontype');
        xy = get(gca,'CurrentPoint');
        kx = xy(1,1);  
        ky = xy(1,2);         
        XLIM = get(Ax_visual, 'xlim');        
        YLIM = get(Ax_visual, 'ylim');
        if queryFlag==1
            disp('Mode 1: disabled clicking!')
            uiresume(f);
        else
            if  gca == Ax_visual  && (kx>=XLIM(1) && kx<=XLIM(2) && ky>=YLIM(1) && ky<=YLIM(2)) 
                Vxy = Xvisual;
                D = (Vxy(:, 1) - kx).^2 + (Vxy(:, 2) - ky).^2;
                [~, jjj] = min(D);  
                uiresume(f);
            else
                switch click_type
                    case {'normal'}  
                        if gca == Ax_EEG_bar  
                            ii = min(max(floor(kx/2)+1, 1), length(labels));   
                            ind_pair = [iEpoch, ii];
                            [~, jjj] = ismember(ind_pair, ind_pairs, 'rows');  
                        elseif gca == Ax_EEG  
                            ii = min(max(floor((kx/Fs)/2)+1, 1), length(labels));   
                            ind_pair = [iEpoch, ii];
                            [~, jjj] = ismember(ind_pair, ind_pairs, 'rows');                         
                        elseif gca == Ax_spec_bar ||gca == Ax_spec{1}||gca == Ax_spec{2}||gca == Ax_spec{3}||gca == Ax_spec{4} ||gca == Ax_spec{5}
                            ii = min(max(floor(kx/2)+1, 1), length(labels));  
                            ind_pair = [iEpoch, ii];
                            [~, jj_] = ismember(ind_pair, ind_pairs, 'rows'); % ind in LUT %
                            jjj = jj_;
                        end
                        uiresume(f);
                end
            end
        end
    end

    function [Ynum, Ystr] = fcn_label2number(Ynum, Ystr)
        if isempty(Ynum)
            Ynum = NaN(length(Ystr), 1);
            Ynum(Ystr == 'Seizure') = 1;
            Ynum(Ystr == 'LPD') = 2;
            Ynum(Ystr == 'GPD') = 3;
            Ynum(Ystr == 'LRDA') = 4;
            Ynum(Ystr == 'GRDA') = 5;
            Ynum(Ystr == 'Other') = 6;
        else
            Ystr = cell(length(Ynum), 1);
            Ystr(Ynum==1) = repmat({'Seizure'}, length(find(Ynum==1)), 1);
            Ystr(Ynum==2) = repmat({'LPD'}, length(find(Ynum==2)), 1);
            Ystr(Ynum==3) = repmat({'GPD'}, length(find(Ynum==3)), 1);
            Ystr(Ynum==4) = repmat({'LRDA'}, length(find(Ynum==4)), 1);
            Ystr(Ynum==5) = repmat({'GRDA'}, length(find(Ynum==5)), 1);
            Ystr(Ynum==6) = repmat({'Other'}, length(find(Ynum==6)), 1);
        end
    end
    
    function fcn_winSize(source, event)
        set(winSizeH, 'Enable', 'off');
        drawnow;
        val = source.Value;
        maps = source.String;       
        winSize = str2double(maps{val}(1:end-4));        
        uiresume(f);
    end

    function fcn_montageFlag(source, event)        
        set(Montage, 'Enable', 'off');
        drawnow;       
        val = source.Value;
        maps = source.String;
        montage = maps{val};       
        uiresume(f);
    end

    function fcn_goNextSPcenter(~, ~)
        set(nextSPcenter, 'Enable', 'off');
        idxReal = [idxReal;{jjj, LUTtopN_current{jj, 5}, datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}];
        idx1_ = idxCPs(find(idxCPs<=jjj, 1, 'last'));  
        idx2_ = idxCPs(find(idxCPs>jjj, 1, 'first'));  
        if isempty(idx2_)  
            tmp_c = jjj;
        else
            tmp_c = floor((idx1_+idx2_)/2); 
        end        
        jj = find(indTopN == tmp_c); 
        LUTtopN_current{jj, 6} = LUTtopN_current{jj, 5};
        LUTtopN{jj, 6} = LUTtopN{jj, 5};
        doneFlag(jj) = 1;
        if queryFlag == 1
            fcn_nextCluster
        else
            fcn_nextSPcenter;
        end
    end

    function Y_bowcpdSmooth = fcn_labelBoWCPDsmooth(~, ~)
        % 2-level smoothing
        % Level 1 - from labeled SP center to all its neighboring+unseen SP centers < R_bow in K^2 distance
        % Level 2 - from SP centers to all its SP memebers            
        idx_labelled = find(doneFlag==1);  
        
        % Step1. Update new labels 
        Y_bowcpdSmooth =  LUT_previous(:, 6);   
        Y_bowcpdSmooth(indTopN(idx_labelled)) = LUTtopN_current(idx_labelled, 6);
        
        % Step2. BoW spread 
        for k = 1:length(idx_labelled)
            kk = idx_labelled(k); 
            
            idx_cc = idxCPCs(kk);   
            ddd_cc = D_bow(kk, :); 

            idx_mm = idxCPCs(setdiff(find(ddd_cc<R_bow), idx_labelled));  
            Y_bowcpdSmooth(idx_mm) = repmat(Y_bowcpdSmooth(idx_cc), length(idx_mm), 1);
        end
 
        % Step3.CPD smooting 
        for k = 1:length(idxCPCs)
            idx1_ = idxCPs(find(idxCPs<=idxCPCs(k), 1, 'last')); 
            idx2_ = -1+idxCPs(find(idxCPs>idxCPCs(k), 1, 'first'));  
            Y_bowcpdSmooth(idx1_:idx2_) = repmat(Y_bowcpdSmooth(idxCPCs(k)), (idx2_-idx1_+1), 1);
        end
        idx_epochEnds = find(isCPs==-1); 
        Y_bowcpdSmooth(idx_epochEnds) = Y_bowcpdSmooth(idx_epochEnds-1);  
    end
 
    function fcn_done(~, ~)      
        idxSeen  = unique([indSeen_; indTopN(doneFlag==1)]);         
        LUT = LUT_previous(:, 1:5);
        LUT(:, 5) = fcn_labelBoWCPDsmooth;       
        timeStamp  = datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF');
        newLUTname = [targetDir,'RaterScores_final.mat'];
        idxReal_idx = cell2mat(idxReal(:, 1));
        [~, idx] = unique(idxReal_idx, 'last');  
        idxReal = idxReal(idx,:);         
        save(newLUTname, 'LUT', 'Y_model', 'idxSeen', 'timeStamp', 'rater', 'idxReal', 'thr_cp', 'idxCPCs');
        choice = questdlg('All set!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        uiresume(f);
                end
    end

    function [icp_, P_, flags1, flags2] = fcn_cpd(sdata, nn)  
        if size(sdata{1, 2}, 2)~=nn
            keyboard
        end
        P = mean(pow2db(cell2mat(sdata(:, 2))+eps), 1); 
        P_ = smooth(P, 10,'sgolay')';
        P_(P_>25) = 25; P_(P_<-10) = -10;
        [icp, ~] = findchangepts(P_, 'Statistic', 'mean','MinThreshold',thr_cp*var(P_));
        flags1 = zeros(nn, 1);
        flags1(icp) = 1;
        flags1(1) = 1;  
        flags1(end) = -1; 
        icp_ = unique([icp, 1, nn]);   
        icp_center = floor((icp_(1:end-1)+icp_(2:end))/2);
        flags2 = zeros(nn, 1);
        flags2(icp_center) = 1;
    end

    function x_final = fcn_parseData(x, itc, i1, i2,  m, n, N)
        if size(x, 2) == n      
            x_final = x;
        else 
            dif = n-size(x, 2);
            if i1 == 1 && i2 ~= N
                x_final = [NaN(20, dif), x];
            elseif i2 == N && i1 ~= 1
                x_final = [x, NaN(m, dif)];
            elseif i2 == N && i1 == 1
                aa = 1- itc-n/2 + 1;
                bb = itc+n/2 - N;
                x_final = [NaN(m, aa), x, NaN(m, bb)];
            end
        end
    end

    function bow = fcn_getBoW(Sdata, NN, K_bow)
        SS = [];
        for i = 1:size(Sdata)
            sdata = Sdata{i, 3}{1};
            nn = length(find(NN==i));           
            if size(sdata{1, 2}, 2)~=nn
                keyboard
            end            
            s = cell2mat(sdata(:, 2))';            
            s = pow2db(s+eps);
            s(s<-10) = -10;
            s(s>25) = 25;
            SS = [SS; s];
        end  
        rng('default')
        bow = kmeans(SS, K_bow);
    end

    function fcn_tm(~, ~)
        set(tm, 'Enable', 'off'); 
        tmp = get(bg,'SelectedObject');
        p = tmp.String;
        disp(p)

        if strcmp(p, 'None')
            tmFlag = -1;
            choice = questdlg('Please pick a sample and label it 1st!','Warning','Ok','Ok');
            switch choice
                case 'Ok'
                    set(0,'CurrentFigure',f)
                    set(f1, 'Visible', 'off');
                    uiresume(f);
            end          
        else
            tmFlag = 1;
            ip = regexp(p,']');
            y_tm = p(ip+2:end);      
            set(f1, 'Visible', 'on');
            set(goBack, 'Enable', 'on');
            for id_r = 1:nP
                set(r1{id_r}, 'Enable', 'on')
            end
            set(ro1, 'Enable', 'on');
            set(bg1,'SelectedObject',ro1)
            for iAx = 1:length(Ax_spec1)
                set(Ax_spec1{iAx}, 'Visible', 'on'); 
            end
            set(Ax_spec_bar1, 'Visible', 'on');   
            set(0,'CurrentFigure',f1)

            indSeen_   = unique([indSeen_; indTopN(doneFlag==1)]); 
            indUnseen_ = setdiff(idxCPCs, indSeen_);       
            if isempty(indUnseen_)
                tmFlag = -1;
                choice = questdlg('All CP center got labeled!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        set(0,'CurrentFigure',f)
                        set(f1, 'Visible', 'off');
                        uiresume(f);
                end 
            else
                idx_tmp = find(ismember(idxCPCs, indUnseen_)==1);
                bow_vec_unseen = bow_vec(idx_tmp, :);  
                idx1_ = idxCPs(find(idxCPs<=jjj, 1, 'last')); 
                idx2_ = idxCPs(find(idxCPs>jjj, 1, 'first')); 
                if isempty(idx2_)  
                    tmp_ = jjj;
                else
                    tmp_ = floor((idx1_+idx2_)/2); 
                end
                bow_vec_current = bow_vec(idxCPCs == tmp_, :);            
                D_= fcn_distChiSq(bow_vec_unseen, bow_vec_current);
                [~, idx_] = min(D_);    
             
                jjj1 = indUnseen_(idx_);  
                iEpoch1 = LUT_previous{jjj1, 1};
                ii1     = LUT_previous{jjj1, 2};
                data1   = DATA{iEpoch1, 2};
                Sdata1  = DATA{iEpoch1, 3}{1}(:, 2);
                stimes1 = DATA{iEpoch1, 3}{2}; 
                sfreqs1 = DATA{iEpoch1, 3}{3};
                totalPower1 = DATA{iEpoch1, 4}{2}; 
                cpds1 = DATA{iEpoch1, 4}{1};            
                idx_local1 = find(cell2mat(LUT_previous(:, 1)) == iEpoch1);             
                labels1  =  LUT_previous(idx_local1, 6);             
                idx_current_local1 = find((doneFlag==1) & (cell2mat(LUTtopN_current(:, 1)) == iEpoch1)); 
                ii_done1 =  cell2mat(LUTtopN_current(idx_current_local1, 2)); 
                labels1(ii_done1) = LUTtopN_current(idx_current_local1, 6);                
                isCPs_local1  = cell2mat(LUT_previous(idx_local1, 3));              
                isCPCs_local1 = cell2mat(LUT_previous(idx_local1, 4));  
                idxCPCs_local1 = find(isCPCs_local1==1);                 
                idxCPs_local1  = find(isCPs_local1~=0);               
                for kk = 1:length(idxCPCs_local1)
                    idx1 = idxCPs_local1(find(idxCPs_local1<=idxCPCs_local1(kk), 1, 'last'));  
                    idx2 = -1+idxCPs_local1(find(idxCPs_local1>idxCPCs_local1(kk), 1, 'first'));  
                    y_cpc = labels1{idxCPCs_local1(kk)};
                    labels1(idx1:idx2) = repmat({y_cpc}, (idx2-idx1+1), 1);
                end
                labels1(end) = labels1(end-1);
                barColors1 = repmat([.8 .8 .8], length(labels1), 1) ;
                Cs_ = flipud(jet(7));
                identifierList_ = categorical(labels1); 
                for kk_ = 1:length(patterns)
                    ind_ = find(identifierList_ == patterns{kk_});
                    barColors1(ind_, :) = repmat(Cs_(kk_,:), length(ind_), 1);
                end       
                tc1 = 2*ii1-1;
                itc1 = tc1*Fs+1;    
                iStart1 = max(1, itc1-Fs*w/2 + 1);            
                iEnd1   = min(size(data1, 2), itc1+Fs*w/2);   
                seg_ = data1(:, iStart1:iEnd1);               
                if size(seg_, 2) == w*Fs  
                    seg1 = seg_;
                else 
                    d = w*Fs-size(seg_, 2);
                    if iStart1 == 1 && iEnd1 ~= size(data1, 2)
                        seg1 = [NaN(20, d), seg_];
                    elseif iEnd1 == size(data1, 2) && iStart1 ~= 1
                        seg1 = [seg_, NaN(20, d)];
                    elseif iEnd1 == size(data1, 2) && iStart1 == 1
                        a = 1- itc1-Fs*w/2 + 1;
                        b = itc1+Fs*w/2 - size(data1, 2);
                        seg1 = [NaN(20, a), seg_, NaN(20, b)];
                    end
                end           
                barColors_eeg_ = barColors1(max(ii1-(w/2-1)/2, 1) : min(length(labels1), ii1+(w/2-1)/2), :);
                barColors_eeg1 = ones(w/2, 3);
                if ii1-(w/2-1)/2<1
                    barColors_eeg1(w/2-size(barColors_eeg_, 1)+1 : w/2, :) = barColors_eeg_;
                elseif ii1+(w/2-1)/2> length(labels1)
                    barColors_eeg1(1: size(barColors_eeg_, 1), :) = barColors_eeg_;
                else
                    barColors_eeg1 = barColors_eeg_;
                end      
                iStart_f1 = max(1, ii1- W/(2*spec_step));   
                iEnd_f1   = min(length(stimes1), ii1+ W/(2*spec_step)-1);
                Sdata_ = cell(4, 1);
                for iii = 1:4
                    S = Sdata1{iii};      
                    seg_f = S(:, iStart_f1:iEnd_f1);
                    if iii == 4
                        barColors_spec_ = barColors1(iStart_f1:iEnd_f1, :);
                        total_power_ = totalPower1(iStart_f1:iEnd_f1);
                    end
                    if size(seg_f, 2) == (W/spec_step)
                        spec1 = seg_f;
                        if iii == 4
                            barColors_spec1 = barColors_spec_;
                            total_power1 = total_power_;
                        end
                    else  
                        d = (W/spec_step) - size(seg_f, 2);
                        if iStart_f1 == 1 && iEnd_f1 ~= length(stimes1)
                            spec1 = [eps+zeros(size(seg_f, 1), d), seg_f];
                            if iii == 4
                                barColors_spec1 = [ones(d, 3); barColors_spec_];
                                total_power1 = [ones(1, d), total_power_];
                            end
                        elseif iEnd_f1 == length(stimes1) && iStart_f1 ~= 1 
                            spec1 = [seg_f, eps+zeros(size(seg_f, 1), d)];
                            if iii == 4
                                barColors_spec1 = [barColors_spec_; ones(d, 3)];
                                total_power1 = [total_power_, ones(1, d)];
                            end
                        elseif iEnd_f1 == length(stimes1) && iStart_f1 == 1 
                            a = 1-(ii1-W/(2*spec_step));
                            b = ii1+ W/(2*spec_step) - length(stimes1)-1;
                            spec1 = [eps+zeros(size(seg_f, 1), a), seg_f, eps+zeros(size(seg_f, 1), b)];
                            if iii == 4
                                barColors_spec1 = [ones(a, 3); barColors_spec_; ones(b, 3)];
                                total_power1 = [ones(1, a), total_power_, ones(1, b)];
                            end
                        end
                    end
                    Sdata_{iii, 1} = spec1;
                end
                Sdata1 = Sdata_;

                timeStampo1 = DATA{iEpoch1, 1}{2};   
                eeg1 = seg1(1:19,:);
                ekg1 = seg1(20, :);               
                tto1 = tc1*Fs - (w/2)*Fs+1;   
                tt11 = tc1*Fs + (w/2)*Fs;    
                tt_1 = tto1:tt11;              
                gap = NaN(1, size(eeg1, 2));
                switch montage
                    case 'L-Bipolar'
                        seg = fcn_bipolar(eeg1);
                        seg_disp1 = [seg(1:4,:); gap; seg(5:8,:); gap; seg(9:12,:); gap; seg(13:16,:); gap; seg(17:18,:); gap; ekg1];
                        channel_withspace1 = channel_withspace_bipolar;
                    case 'Average' 
                        seg = eeg - repmat(mean(eeg, 1), size(eeg, 1), 1);
                        seg_disp1 = [seg(1:8,:); gap; seg(9:11,:); gap; seg(12:19,:); gap; ekg1];
                        channel_withspace1 = channel_withspace_average;
                    case 'Monopolar'
                        seg =  eeg;
                        seg_disp1 = [seg(1:8,:); gap; seg(9:11,:); gap; seg(12:19,:); gap; ekg1];
                        channel_withspace1 = channel_withspace_monopolar;
                end
                M1 = size(seg_disp1, 1);
                DCoff1 = repmat(flipud((1:M1)'), 1, size(seg_disp1, 2));
                timeStamps1 = datestr(timeStampo1 + seconds(round(tto1/Fs:2:tt11/Fs)), 'hh:MM:ss');
                set(f1,'CurrentAxes',Ax_EEG_bar1);cla(Ax_EEG_bar1)
                colors_ = reshape(barColors_eeg1 , 1, size(barColors_eeg1, 1), 3); 
                image(Ax_EEG_bar1, round(tt_1(1)/Fs)+(1:2:14), 1, colors_);
                set(Ax_EEG_bar1,'ytick',[],'xtick',[])
                box(Ax_EEG_bar1, 'on'); 
                set(f1,'CurrentAxes',Ax_EEG1);cla(Ax_EEG1)
                hold(Ax_EEG1, 'on')                
                    for iSec = 1:round((tt11-tto1+1)/Fs)
                        ta = tto1 + Fs*(iSec-1);
                        if iSec== 7 || iSec==9
                            line(Ax_EEG1,[ta ta],  [0 M1+1], 'linestyle', '--', 'color', 'r')
                        else
                            line(Ax_EEG1,[ta ta],  [0 M1+1], 'linestyle', '--', 'color', [.5 .5 .5])
                        end
                    end
                    plot(Ax_EEG1, tt_1, zScale*seg_disp1(1:end-1,:)+DCoff1(1:end-1,:),'k');
                    ekg_ = seg_disp1(end,:);
                    ekg_1 = (ekg_-nanmean(ekg_))/(eps+nanstd(ekg_));
                    plot(Ax_EEG1, tt_1, .2*ekg_1+DCoff1(end,:),'m');
                    box(Ax_EEG1, 'on');
                    axis(Ax_EEG1, 'on');
                    set(Ax_EEG1, 'ytick',1:M1,'yticklabel',fliplr(channel_withspace1),'box','on', 'ylim', [0 M1+1], 'xlim', [tto1 tt11], 'xtick',round(tt_1(1):2*Fs:tt_1(end)),'xticklabel',timeStamps1)
 
                    dt = tt11-tto1+1;
                    a = round(dt*4/5);
                    xa1 = tto1+[a a+Fs-1];        
                    ya1 = [3 3];
                    xa2 = tto1+[a a];        
                    ya2 = ya1+[0 100*zScale];
                    text(Ax_EEG1, xa1(1)-.7*a/10,    mean(ya2), '100\muV','Color', 'r','FontSize',12);
                    text(Ax_EEG1, mean(xa1)-0.3*a/10, 2.5, '1 sec','Color', 'r','FontSize',12);        
                    line(Ax_EEG1, xa1,ya1, 'LineWidth', 2, 'Color','r');
                    line(Ax_EEG1, xa2,ya2, 'LineWidth', 2, 'Color','r');
                hold(Ax_EEG1, 'off')  
               
                to1 = (tc1-1)-(winSize*60/2)+1; 
                t11 = (tc1-1)+(winSize*60/2);              
                S_x1 = to1:stepSize:t11;   
                S_y1 = sfreqs1;
                for iReg = 1:4                  
                    set(f1,'CurrentAxes',Ax_spec1{iReg}); cla(Ax_spec1{iReg})
                    spec = Sdata1{iReg};colormap jet;
                    imagesc(Ax_spec1{iReg}, S_x1, S_y1, pow2db(spec), col);  axis(Ax_spec1{iReg},'xy'); 
                    box(Ax_spec1{iReg}, 'on');  
                    set(Ax_spec1{iReg},'xtick',[])   
                    hold(Ax_spec1{iReg}, 'on')
                        plot(Ax_spec1{iReg}, [tc1  tc1], [S_y1(1) S_y1(end)], 'k--','linewidth',1);
                    hold(Ax_spec1{iReg}, 'off')                   
                    xx = get(Ax_spec1{iReg}, 'yticklabel'); xx{end} = spatialRegs{iReg};
                    xlim([tc1-winSize/2*60+1, tc1+winSize/2*60])
                    set(Ax_spec1{iReg}, 'yticklabel',xx)
                    ylabel(Ax_spec1{iReg}, 'Freq (Hz)')
                end
                
                set(f1,'CurrentAxes',Ax_spec1{5});     cla(Ax_spec1{5})
                hold(Ax_spec1{5}, 'on')  
                    if length(S_x1)~=length(total_power1)
                        keyboard
                    end
                    plot(Ax_spec1{5}, S_x1, total_power1, 'g-', 'linewidth', 1)
                    cpds_tmp_ = stimes1(cpds1)-1;
                    idx_tmp = find(cpds_tmp_>=tc1-winSize/2*60+1 & cpds_tmp_<=tc1+winSize/2*60);
                    cpds_tmp1 = cpds_tmp_(max(1, min(idx_tmp)-1) : min(length(cpds_tmp_), max(idx_tmp)+1));
                    for icpd = 1:length(cpds_tmp1)-1
                        aa_ = cpds_tmp1(icpd);
                        bb_ = cpds_tmp1(icpd+1);
                        plot(Ax_spec1{5}, [aa_, aa_], [-10 25], 'm-', 'linewidth', .5)
                        cpd_mean = mean(totalPower1((aa_/2+1):(bb_/2+1)));
                        plot(Ax_spec1{5}, [aa_, bb_], [cpd_mean, cpd_mean], 'b-', 'linewidth', .5)
                    end
                    plot(Ax_spec1{5}, [tc1  tc1], [-10 25], 'k--','linewidth',1);                    
                    tt_1 = timeStampo1 + seconds(to1:10*60:t11);
                    set(Ax_spec1{5}, 'xtick', to1:10*60:t11, 'xticklabel', datestr(tt_1, 'hh:MM:ss'),...
                        'xlim',[tc1-winSize/2*60+1, tc1+winSize/2*60],...
                        'ylim', [-10 25],...
                        'ytick', [], 'yticklabel', []);
                    box(Ax_spec1{5}, 'on')
                hold(Ax_spec1{5}, 'off')  
                
                set(f1,'CurrentAxes',Ax_spec1{end});
                cla(Ax_spec1{end})
                plot(Ax_spec1{end}, tc1, 0,'rv','markersize', 8,'MarkerFaceColor','r'); 
                set(Ax_spec1{end},'xtick',[],'xlim',get(Ax_spec1{end-1},'xlim'), 'ylim', [0 .5]) 
                axis(Ax_spec1{end}, 'off');
                
                set(f1,'CurrentAxes',Ax_spec_bar1);cla(Ax_spec_bar1)
                colors_ = reshape(barColors_spec1 , 1, size(barColors_spec1, 1), 3); 
                image(Ax_spec_bar1, S_x1, 1, colors_);
                set(Ax_spec_bar1,'xtick',[],'ytick',[])
                xlim([tc1-winSize/2*60+1, tc1+winSize/2*60])
                box(Ax_spec_bar1, 'on')
            end
        end      
        uiresume(f);
    end

    function fcn_tmdecision(source, event)
        for m = 1:nP
            set(r1{m}, 'Enable', 'off');
        end
        set(ro1, 'Enable', 'off');
        drawnow;      
        tmFlag = 1; set(f1, 'Visible', 'on');     
        p = event.NewValue.String;
        ip = regexp(p,']');
        p = p(ip+2:end);
        idxReal = [idxReal;{jjj1, p, datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}];       
        idx1_ = idxCPs(find(idxCPs<=jjj1, 1, 'last'));  
        idx2_ = idxCPs(find(idxCPs>jjj1, 1, 'first'));  
        if isempty(idx2_)  
            tmp_c = jjj1;
        else
            tmp_c = floor((idx1_+idx2_)/2); 
        end
        if isempty(find(indTopN == tmp_c, 1))         
            keyboard
        else
            jj_c = find(indTopN == tmp_c); 
        end        
        LUTtopN_current{jj_c, 6} = p;
        doneFlag(jj_c) = 1;       
        if strcmp(y_tm, p)
            fcn_tm();
        else 
            tmFlag = 0; set(f1, 'Visible', 'off');
            uiresume(f);
        end 
    end

    function fcn_tmdecisionkeys(ind_)       
        p = patterns{ind_};
        idxReal = [idxReal;{jjj1, p, datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}];
        idx1_ = idxCPs(find(idxCPs<=jjj1, 1, 'last'));  
        idx2_ = idxCPs(find(idxCPs>jjj1, 1, 'first'));  
        if isempty(idx2_) 
            tmp_c = jjj1;
        else
            tmp_c = floor((idx1_+idx2_)/2); 
        end     
        if isempty(find(indTopN == tmp_c, 1))          
            keyboard;
        else
            jj_c = find(indTopN == tmp_c);  
        end       
        LUTtopN_current{jj_c, 6} = p;
        doneFlag(jj_c) = 1;
        if strcmp(y_tm, p)
            fcn_tm();
        else
            tmFlag = 0; set(f1, 'Visible', 'off');
            uiresume(f);
        end
    end

    function fcn_goBack(varargin)
        set(goBack, 'Enable', 'off');
        tmFlag = 0;
        set(f1, 'Visible', 'off');
        uiresume(f);
    end

    function D = fcn_distChiSq(X, Y)
        % The chi-squared distance between two vectors is defined as:
        % d(x,y) = sum( (xi-yi)^2 / (xi+yi) ) / 2;
        % The chi-squared distance is useful when comparing histograms.
        mm = size(X,1);  
        nn = size(Y,1);
        mOnes = ones(1,mm); D = zeros(mm,nn);
        for i=1:nn
            yi = Y(i,:);  
            yiRep = yi(mOnes, :);
            s = yiRep + X;    
            ddd = yiRep - X;
            D(:,i) = sum(ddd.^2./(s+eps), 2);
        end
        D = D/2;
    end

    function fcn_nextSPcenter(varargin)       
        jj = find(isnan(doneFlag), 1);
        if isempty(jj)
            jj = length(indTopN);
            choice = questdlg('Reach the last cluster!','Warning','Ok','Ok');
            switch choice
                case 'Ok'
                    jjj = size(LUT_previous, 1);
            end 
        else
            jjj = indTopN(jj);           
        end
        uiresume(f);    
    end

    function fcn_radiusSize(source, event)
        set(radiusSizeH, 'Enable', 'off');
        drawnow;
        val = source.Value;
        maps = source.String;        
        R_bow = str2double(maps{val});       
        uiresume(f);
    end

    function fcn_acceptNnext(~, ~)
        set(acceptNnext, 'Enable', 'off');       
        drawnow;   
        tmFlag = 0; set(f1, 'Visible', 'off');
        jjj_1 = jjj-(winSize*60/2)+1;
        jjj_2 = jjj+(winSize*60/2);
        idx_currentWin = find(idxCPCs>=jjj_1 & idxCPCs<=jjj_2);
        doneFlag(idx_currentWin) = 1;
        LUTtopN_current(idx_currentWin, 6) = LUTtopN_current(idx_currentWin, 5);
        idxReal = [idxReal; [num2cell(idxCPCs(idx_currentWin)), LUTtopN_current(idx_currentWin, 5), repmat({datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}, length(idx_currentWin), 1)]];
        jjj_ = jjj+(winSize*60/2);
        if jjj_>size(LUT_previous, 1)
            choice = questdlg('Reach the end!','Warning','Ok','Ok');
            switch choice
                case 'Ok'
                    jjj = size(LUT_previous, 1);
            end
        else
            jjj = jjj_;
        end
        uiresume(f);
    end

    function [idx_mediods, idx_member2medoid] = fcn_kmediodDistChiSq(X, kMediods)

        [idx, C] = kmedoids(X, kMediods, 'Distance', @fcn_distChiSq);       
        idx_mediods = NaN(kMediods, 1);
        idx_member2medoid = NaN(length(idx), 1);
        for k=1:kMediods
            kk = find(ismember(X, C(k, :), 'rows'));
            idx_mediods(k) = kk(1);
            idx_member2medoid(idx==k) =  kk(1);
        end      
        idx_mediods = sort(idx_mediods);
    end

    function Y_cpdSmooth = fcn_labelClusterCPDsmooth(~, ~) 
        % Step1. Update new labels 
        Y_cpdSmooth = LUT_previous(:, 6);%
        Y_cpdSmooth(indTopN(doneFlag==1)) = LUTtopN_current(doneFlag==1, 6);
        
        % Step2. Cluster smooth 
        for k = 1:length(idx_mediods_cpc)
            idx_cc = idxCPCs(idx_mediods_cpc(k));
            idx_cm = idxCPCs(idx_member2medoid_cpc==idx_mediods_cpc(k));          
            Y_cpdSmooth(idx_cm) = repmat(Y_cpdSmooth(idx_cc), length(idx_cm), 1);
        end
        
        % Step3.CPD smooting 
        for k = 1:length(idxCPCs)
            idx1_ = idxCPs(find(idxCPs<=idxCPCs(k), 1, 'last'));  
            idx2_ = -1+idxCPs(find(idxCPs>idxCPCs(k), 1, 'first'));  
            Y_cpdSmooth(idx1_:idx2_) = repmat(Y_cpdSmooth(idxCPCs(k)), (idx2_-idx1_+1), 1);
        end
        idx_epochEnds = find(isCPs==-1); 
        Y_cpdSmooth(idx_epochEnds) = Y_cpdSmooth(idx_epochEnds-1);  
    end

    function [cc, y] = fcn_colorClusterCPDsmooth(~, ~)
        % Step1. Update new labels 
        y = LUT_previous(:, 6);%  
        y(indTopN(doneFlag==1)) = LUTtopN_current(doneFlag==1, 6);
        y(indTopN(isnan(doneFlag))) = repmat({'Unknown'}, sum(isnan(doneFlag)), 1);
        
        % Step2. Cluster smooth 
        for k = 1:length(idx_mediods_cpc)
            idx_cc = idxCPCs(idx_mediods_cpc(k));
            idx_cm = idxCPCs(idx_member2medoid_cpc==idx_mediods_cpc(k));           
            y(idx_cm) = repmat(y(idx_cc), length(idx_cm), 1);
        end
        
        % Step3.CPD smooth
        for k = 1:length(idxCPCs)
            idx1_ = idxCPs(find(idxCPs<=idxCPCs(k), 1, 'last'));  
            idx2_ = -1+idxCPs(find(idxCPs>idxCPCs(k), 1, 'first'));  
            y(idx1_:idx2_) = repmat(y(idxCPCs(k)), (idx2_-idx1_+1), 1);
        end
        idx_epochEnds = find(isCPs==-1); 
        y(idx_epochEnds) = y(idx_epochEnds-1);  
 
        % Step4.turn into colors %
        cc = repmat([.8 .8 .8], size(LUT_previous, 1), 1);
        Cs = flipud(jet(7));
        scr = categorical(y);          
        for k = 1:length(patterns)
            idx = find(scr == patterns{k});
            cc(idx, :) = repmat(Cs(k,:), length(idx), 1);
        end
    end

    function fcn_nextCluster(varargin) 
        jj = find(isnan(doneFlag), 1);
        if isempty(jj)  
            queryFlag = 2;
            set(done, 'Visible', 'on');
            set(nextSPcenter, 'Visible', 'on');
            set(acceptNnext, 'Visible', 'on');
            set(tm, 'Visible', 'on');
            set(allSZ, 'Visible', 'on'); set(allLPD, 'Visible', 'on'); set(allGPD, 'Visible', 'on'); set(allLRDA, 'Visible', 'on'); set(allGRDA, 'Visible', 'on'); set(allOther, 'Visible', 'on');
            idxSeen = indTopN; 
            LUT_previous(:, 5) = fcn_labelClusterCPDsmooth;
            LUT_previous(:, 6) = LUT_previous(:, 5);            
            LUT = LUT_previous(:, 1:5); 
            timeStamp  = datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF');
            newLUTname = [targetDir,'RaterScores_clusters.mat'];
            idxReal_idx = cell2mat(idxReal(:, 1));
            [~, idx] = unique(idxReal_idx, 'last'); 
            idxReal = idxReal(idx,:);           
            save(newLUTname, 'LUT', 'Y_model', 'idxSeen', 'timeStamp', 'rater', 'idxReal', 'thr_cp', 'idxCPCs', 'idx_mediods_cpc', 'idx_member2medoid_cpc');
  
            indTopN = idxCPCs;
            LUTtopN = LUT_previous(indTopN, :);            
            N = length(indTopN);                              
            LUTtopN_current = cell(N, 6);                       
            LUTtopN_current(:, 1:5) = LUTtopN(:, [1:4, 6]);    
            doneFlag = NaN(N, 1);
            doneFlag(ismember(indTopN, idxSeen)) = 1; 
            LUTtopN_current(doneFlag==1, 6) = LUTtopN_current(doneFlag==1, 5);           
            jj = find(isnan(doneFlag), 1);            
            if ~isempty(jj)
                choice = questdlg('Enter sequential mode!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        jjj = indTopN(jj);
                end 
            else
                choice = questdlg('All samples got labelled!','Warning','Ok','Ok');
                switch choice
                    case 'Ok'
                        fcn_done;                         
                end 
            end
        else
            jjj = indTopN(jj);         
        end
        uiresume(f);    
    end

    function fcn_allSZ(~, ~)
        set(allSZ, 'Enable', 'off');
        drawnow; 
        iEpoch_local = cell2mat(LUT_previous(:, 1));         
        iStart_local = max(1, jjj-60/2*winSize/2 + 1);        
        iEnd_local   = min(nSamples, jjj+60/2*winSize/2);            
        idx1_ = idxCPs(find(idxCPs<=iStart_local, 1, 'last'));  
        idx2_ = idxCPs(find(idxCPs>iStart_local, 1, 'first'));  
        if isempty(idx2_)  
            tmp_1 = iStart_local;
        else
            tmp_1 = floor((idx1_+idx2_)/2); 
        end        
        idx1_ = idxCPs(find(idxCPs<=iEnd_local, 1, 'last'));  
        idx2_ = idxCPs(find(idxCPs>iEnd_local, 1, 'first'));  
        if isempty(idx2_)  
            tmp_2 = iEnd_local;
        else
            tmp_2 = floor((idx1_+idx2_)/2); 
        end        
        idx = find(isCPCs==1 & iEpoch_local==iEpoch);  
        idx = idx(idx>=tmp_1 & idx<=tmp_2);
        idxReal = [idxReal; [num2cell(idx), repmat({'Seizure', datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}, length(idx), 1)]];        
        for i_ = 1:length(idx)
            tmp_ = idx(i_);
            if isempty(find(indTopN == tmp_, 1))  
                keyboard  
            else
                jj_ = find(indTopN == tmp_);     
            end            
            LUTtopN_current{jj_, 6} = 'Seizure';
            doneFlag(jj_) = 1;
        end
        fcn_nextSPcenter;
    end

    function fcn_allLPD(~, ~)
        set(allLPD, 'Enable', 'off');
        drawnow;
        iEpoch_local = cell2mat(LUT_previous(:, 1));         
        iStart_local = max(1, jjj-60/2*winSize/2 + 1);       
        iEnd_local   = min(nSamples, jjj+60/2*winSize/2);            
        idx1_ = idxCPs(find(idxCPs<=iStart_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iStart_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_1 = iStart_local;
        else
            tmp_1 = floor((idx1_+idx2_)/2); 
        end      
        idx1_ = idxCPs(find(idxCPs<=iEnd_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iEnd_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_2 = iEnd_local;
        else
            tmp_2 = floor((idx1_+idx2_)/2); 
        end       
        idx = find(isCPCs==1 & iEpoch_local==iEpoch); 
        idx = idx(idx>=tmp_1 & idx<=tmp_2);       
        idxReal = [idxReal; [num2cell(idx), repmat({'LPD', datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}, length(idx), 1)]];
        for i_ = 1:length(idx)
            tmp_ = idx(i_);
            if isempty(find(indTopN == tmp_, 1)) 
                keyboard  
            else
                jj_ = find(indTopN == tmp_);     
            end          
            LUTtopN_current{jj_, 6} = 'LPD';
            doneFlag(jj_) = 1;
        end
        fcn_nextSPcenter;
    end

    function fcn_allGPD(~, ~)
        set(allGPD, 'Enable', 'off');
        drawnow;
        iEpoch_local = cell2mat(LUT_previous(:, 1));      
        iStart_local = max(1, jjj-60/2*winSize/2 + 1);       
        iEnd_local   = min(nSamples, jjj+60/2*winSize/2);           
        idx1_ = idxCPs(find(idxCPs<=iStart_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iStart_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_1 = iStart_local;
        else
            tmp_1 = floor((idx1_+idx2_)/2); 
        end       
        idx1_ = idxCPs(find(idxCPs<=iEnd_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iEnd_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_2 = iEnd_local;
        else
            tmp_2 = floor((idx1_+idx2_)/2); 
        end       
        idx = find(isCPCs==1 & iEpoch_local==iEpoch); 
        idx = idx(idx>=tmp_1 & idx<=tmp_2);
        idxReal = [idxReal; [num2cell(idx), repmat({'GPD', datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}, length(idx), 1)]];       
        for i_ = 1:length(idx)
            tmp_ = idx(i_);
            if isempty(find(indTopN == tmp_, 1)) 
                keyboard  
            else
                jj_ = find(indTopN == tmp_);     
            end            
            LUTtopN_current{jj_, 6} = 'GPD';
            doneFlag(jj_) = 1;
        end
        fcn_nextSPcenter;
    end

    function fcn_allLRDA(~, ~)
        set(allLRDA, 'Enable', 'off');
        drawnow;
        iEpoch_local = cell2mat(LUT_previous(:, 1));        
        iStart_local = max(1, jjj-60/2*winSize/2 + 1);       
        iEnd_local   = min(nSamples, jjj+60/2*winSize/2);           
        idx1_ = idxCPs(find(idxCPs<=iStart_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iStart_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_1 = iStart_local;
        else
            tmp_1 = floor((idx1_+idx2_)/2); 
        end        
        idx1_ = idxCPs(find(idxCPs<=iEnd_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iEnd_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_2 = iEnd_local;
        else
            tmp_2 = floor((idx1_+idx2_)/2); 
        end       
        idx = find(isCPCs==1 & iEpoch_local==iEpoch); 
        idx = idx(idx>=tmp_1 & idx<=tmp_2); 
        idxReal = [idxReal; [num2cell(idx), repmat({'LRDA', datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}, length(idx), 1)]];       
        for i_ = 1:length(idx)
            tmp_ = idx(i_);
            if isempty(find(indTopN == tmp_, 1)) 
                keyboard  
            else
                jj_ = find(indTopN == tmp_);     
            end            
            LUTtopN_current{jj_, 6} = 'LRDA';
            doneFlag(jj_) = 1;
        end
        fcn_nextSPcenter;
    end

    function fcn_allGRDA(~, ~)
        set(allGRDA, 'Enable', 'off');
        drawnow;      
        iEpoch_local = cell2mat(LUT_previous(:, 1));         
        iStart_local = max(1, jjj-60/2*winSize/2 + 1);       
        iEnd_local   = min(nSamples, jjj+60/2*winSize/2);            
        idx1_ = idxCPs(find(idxCPs<=iStart_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iStart_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_1 = iStart_local;
        else
            tmp_1 = floor((idx1_+idx2_)/2); 
        end       
        idx1_ = idxCPs(find(idxCPs<=iEnd_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iEnd_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_2 = iEnd_local;
        else
            tmp_2 = floor((idx1_+idx2_)/2); 
        end
        idx = find(isCPCs==1 & iEpoch_local==iEpoch); 
        idx = idx(idx>=tmp_1 & idx<=tmp_2);
        idxReal = [idxReal; [num2cell(idx), repmat({'GRDA', datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}, length(idx), 1)]];      
        for i_ = 1:length(idx)
            tmp_ = idx(i_);
            if isempty(find(indTopN == tmp_, 1)) 
                keyboard  
            else
                jj_ = find(indTopN == tmp_);     
            end            
            LUTtopN_current{jj_, 6} = 'GRDA';
            doneFlag(jj_) = 1;
        end
        fcn_nextSPcenter;
    end

    function fcn_allOther(~, ~)
        set(allOther, 'Enable', 'off');
        drawnow;
        iEpoch_local = cell2mat(LUT_previous(:, 1));         
        iStart_local = max(1, jjj-60/2*winSize/2 + 1);       
        iEnd_local   = min(nSamples, jjj+60/2*winSize/2);            
        idx1_ = idxCPs(find(idxCPs<=iStart_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iStart_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_1 = iStart_local;
        else
            tmp_1 = floor((idx1_+idx2_)/2); 
        end       
        idx1_ = idxCPs(find(idxCPs<=iEnd_local, 1, 'last')); 
        idx2_ = idxCPs(find(idxCPs>iEnd_local, 1, 'first')); 
        if isempty(idx2_) 
            tmp_2 = iEnd_local;
        else
            tmp_2 = floor((idx1_+idx2_)/2); 
        end      
        idx = find(isCPCs==1 & iEpoch_local==iEpoch); 
        idx = idx(idx>=tmp_1 & idx<=tmp_2);
        idxReal = [idxReal; [num2cell(idx), repmat({'Other', datestr(now, 'yyyy-mmm-dd HH:MM:ss.FFF')}, length(idx), 1)]];
        for i_ = 1:length(idx)
            tmp_ = idx(i_);
            if isempty(find(indTopN == tmp_, 1)) 
                keyboard  
            else
                jj_ = find(indTopN == tmp_);     
            end
            LUTtopN_current{jj_, 6} = 'Other';
            doneFlag(jj_) = 1;
        end
        fcn_nextSPcenter;
    end
end