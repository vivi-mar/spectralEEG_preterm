% driver to process an EEG dataset

warning off

% ****************************************
% !!! add EEGLAB to the path first
% ****************************************

[FN, PN]  = uigetfile('/*******/*.edf', 'select .EDF file to extract subject DATA', 'multiselect', 'on');  % add the link to your folder
 
 for K = 1 : size(FN,2)
    entirefile =fullfile(PN,FN{K}) 
    fid1 = fopen(entirefile);
    
    [Path, FileName, extension] = fileparts(entirefile);

    EEG         = pop_biosig(strcat(PN, FN{K}), 'importevent', 'off');
    EEG.setname = FN{K};
    EEG         = eeg_checkset(EEG);

        for i = 12:size(EEG.chanlocs, 1)
            EEG.chanlocs(i).labels = EEG.chanlocs(i).labels(5:end);
        end

    EEG         = pop_chanedit(EEG, 'lookup', 'standard-10-5-cap385.elp');

    EEG         = pop_select(EEG, 'channel', {'Fp1';'Fp2';'T3';'C3';'C4';'T4';'O1';'O2'}); % change the labels based on the .edf file

    % Filter EEG data imported from EDF

    sensors = 1:size(EEG.data, 1); % ID of channels to be filtered


    % Notch filter 50Hz EEG

    wo     = 50/(EEG.srate/2);
    bw     = wo/50;
    [b, a] = iirnotch(wo, bw);
    h = waitbar(0, 'EEG Performing IIR notch filter...');
        for i = 1:size(sensors, 2)
            EEG.signal(sensors(i), :) = filtfilt(b, a, double(EEG.data(sensors(i), :)));
            waitbar(i/size(sensors,2))
        end
    close(h)
    clear b a h

    % Bandpass filter EEG

    lp = 32;
    hp = .1;

    w2 = 2*lp/EEG.srate;
    w  = w2;
    [B, A] = butter(4, w, 'low');
    w2 = 2*hp/EEG.srate;
    w  = w2;
    [C, D] = butter(4, w, 'high');

    h = waitbar(0,['EEG Performing ', num2str(hp), '-', num2str(lp), ' Bandpass Filter...']);

        for kk = 1:size(sensors, 2)
            waitbar(kk/size(sensors, 2))
            EEG.data(sensors(kk), :) = filtfilt(C, D, filtfilt(B, A, double(EEG.data(sensors(kk), :))));
        end
    close(h)
    clear w w2 A B C D kk a h

    % Plot EEG % 
    %if you want to see the result of filtering uncomment the following line
    % eegplot(EEG.data, 'srate', EEG.srate, 'winlength', 60)

    
    %% Interpolate bad channels selected visually
    EEG.etc.ignoredtime = 10; % input('Seconds to ignore: '); % do not include first part of recording in batch selection (exclude 10s)

    tempdata = EEG.data;
    tempch   = [];
    distr    = [];
    h = waitbar(0, 'Looking for bad channels...');
    for kk = 1:size(EEG.data, 1)
        waitbar(kk/size(EEG.data, 1))
        if  (max(tempdata(kk, (EEG.etc.ignoredtime*EEG.srate:end))) - min(tempdata(kk, (EEG.etc.ignoredtime*EEG.srate:end)))) > 1000 || ...
             max(tempdata(kk, (EEG.srate:end))) - min(tempdata(kk, (EEG.srate:end))) < 20 % max(tempdata(kk,:))>=150 || min(tempdata(kk,:))<=-150
            tempch = [tempch, kk];
        end
        temp  = max(tempdata(kk, (EEG.etc.ignoredtime*EEG.srate:end))) - min(tempdata(kk, (EEG.etc.ignoredtime*EEG.srate:end)));
        distr = [distr, temp];
    end
    close (h)

    % Find additional badchannels using signal distribution params

    distr(tempch) = NaN;
    otherbad      = find(distr >= nanmean(distr) + 3.5*nanstd(distr));

    disp(['Method 1: ', num2str(tempch),' <---> Method 2: ', num2str(otherbad)]); 

    % Before going ahead manually validate automatic bad channel selection and include manually spotted bad chs

            if size(tempch,2) == EEG.nbchan
                eegplot(EEG.data, 'winlength', 60)
                bad_elec = input('Channels to reject: ');
                EEG.etc.badchannels = bad_elec;
                badchannels         = bad_elec;
                clear bad_elec
            else
                EEG.etc.badchannels = union(tempch, otherbad);
                badchannels         = union(tempch, otherbad);
            end
    clear kk j tempdata h

    method  = 'spherical';
    EEGtemp = eeg_emptyset;

    EEGtemp.specdata    = [];
    EEGtemp.icachansind = [];
    EEGtemp.specicaact  = [];
    EEGtemp.reject      = [];
    EEGtemp.stats       = [];

    EEGtemp.nbchan      = size(EEG.data, 1);
    EEGtemp.data        = EEG.data;
    EEGtemp.chanlocs    = EEG.chanlocs;
    EEGtemp.trials      = size(EEG.data, 3);
    EEGtemp.pnts        = size(EEG.data, 2);
    EEGtemp.srate       = EEG.srate;
    EEGtemp.badchannels = EEG.etc.badchannels;

        if ~isempty(EEGtemp.badchannels)
            EEGOUT = eeg_interp(EEGtemp, EEGtemp.badchannels, method);
            EEGtemp.data = EEGOUT.data;
        end

    EEG.data = EEGtemp.data;
    clear EEGOUT EEGtemp method
    
    Pathfilt = '/Volumes/****/EEG_filtered/'; % !! to be changed based on the selected folder for results
    FileName = strcat(FileName);
    save(strcat(Pathfilt, FileName), 'EEG');
 end
 
 warning on
