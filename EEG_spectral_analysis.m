% driver to process an EEG dataset
clear all 
close all

warning off

% ****************************************
% load data (add EEGLAB to the path first)
% ****************************************

[FN, PN]  = uigetfile('/*******/EEG_filtered/*.mat', 'multiselect', 'on'); % paste the location of the folder containing the filtered files
fid = fopen('total_spectra.csv','w');
for K = 1 : size(FN,2)

      inpfilename1 = FN{K}
      entirefile =fullfile(PN,FN{K})
      [Path, FileName, extension] = fileparts(entirefile)
      load(inpfilename1);

    % *****************
    % spectral analysis
    % *****************

    N_win1     = 4*EEG.srate;
    N_overlap1 = N_win1/2;
    N_fft1     = 4096;
    N_pts1     = size(EEG.data, 2);
    N_ch1      = size(EEG.data, 1);

    
    GR1.gruppi1.FRO=double(EEG.data(1,:)-EEG.data(2,:));
    GR1.gruppi1.CEN=double(EEG.data(4,:)-EEG.data(5,:));
    GR1.gruppi1.OCC=double(EEG.data(7,:)-EEG.data(8,:));
    GR1.gruppi1.CT_DX=double(EEG.data(5,:)-EEG.data(6,:));
    GR1.gruppi1.CT_SX=double(EEG.data(4,:)-EEG.data(3,:));
    GR1.gruppi1.TOT=mean(double(EEG.data(:,:)));

    matrice1 = cell2mat(cellfun(@double,struct2cell(GR1.gruppi1),'UniformOutput',false));
   
  
    %%

    [Pxx1, F1] = pwelch( matrice1', hamming(N_win1), N_overlap1, N_fft1, EEG.srate, 'power');
    DF1 = F1(2) - F1(1);
    
    Pxxfinal1= Pxx1;
    grpn1=size(Pxxfinal1,2);
    grpn_name1= fieldnames(GR1.gruppi1); 
    
    %%
    % *******
    % delta 1 
    % *******

    I_lo_delta_1 = find(F1 == 0.5);
    I_hi_delta_1 = find(F1 == 1);

    % *******
    % delta 2 
    % *******

    I_lo_delta_2 = find(F1 > 1, 1, 'first');
    I_hi_delta_2 = find(F1 == 4);

    % *****
    % theta 
    % *****

    I_lo_theta = find(F1 > 4, 1, 'first');
    I_hi_theta = find(F1 == 8);

    % *****
    % alpha
    % *****

    I_lo_alpha = find(F1 > 8, 1, 'first');
    I_hi_alpha = find(F1 == 13);

    % ****
    % beta
    % ****

    I_lo_beta = find(F1 > 13, 1, 'first');
    I_hi_beta = find(F1 == 30);
    
    % ****
    % total
    % ****

    I_lo_total = find(F1 > 0.5, 1, 'first');
    I_hi_total = find(F1 == 30);
    

    % *******
    % summary
    % *******

    abs_power1 = trapz(Pxxfinal1)*DF1;
    power1     = zeros(N_fft1/2+1, grpn1);
    
    %%
    for i = 1:grpn1
        power1(:, i) = cumtrapz(Pxxfinal1(:, i)).*(DF1/abs_power1(i));
    end

    spect1.abs_D1 = trapz(Pxxfinal1(I_lo_delta_1:I_hi_delta_1, :))*DF1; 
    spect1.abs_D2 = trapz(Pxxfinal1(I_lo_delta_2:I_hi_delta_2, :))*DF1;
    spect1.abs_T  = trapz(Pxxfinal1(I_lo_theta:I_hi_theta, :))*DF1;
    spect1.abs_A  = trapz(Pxxfinal1(I_lo_alpha:I_hi_alpha, :))*DF1;
    spect1.abs_B  = trapz(Pxxfinal1(I_lo_beta:I_hi_beta, :))*DF1;
   
    abs_pwr_mat1 = cell2mat(cellfun(@double,struct2cell(spect1),'UniformOutput',false));
    A1=transpose(abs_pwr_mat1);
    abs_band_name1=fieldnames(spect1);
    
    abs_power_total1 = trapz(Pxxfinal1(I_lo_total:I_hi_total, :))*DF1;
    abs_pwrlog_total1 = 10*log10(trapz(Pxxfinal1(I_lo_total:I_hi_total, :))*DF1);
    
    perc1.rel_D1 = spect1.abs_D1./abs_power_total1.*100;
    perc1.rel_D2 = spect1.abs_D2./abs_power_total1.*100;
    perc1.rel_T  = spect1.abs_T./abs_power_total1.*100;
    perc1.rel_A  = spect1.abs_A./abs_power_total1.*100;
    perc1.rel_B  = spect1.abs_B./abs_power_total1.*100;

    rel_pwr_mat1 = cell2mat(cellfun(@double,struct2cell(perc1),'UniformOutput',false));
    R1=transpose(rel_pwr_mat1);
    rel_band_name1=fieldnames(perc1);

    SEF1 = zeros(grpn1, 1);
    for i = 1:grpn1
        I      = find(power1(:, i) >= 0.95, 1, 'first');
        SEF1(i) = F1(I);
    end
    
    ABSband_name=fieldnames(spect1);
    RELband_name=fieldnames(perc1);



 
%% writes report

         fprintf(fid,';%s',[])
       for b=1:length(grpn_name1);
        for c=1:length(abs_band_name1);
            fprintf(fid,';%s',[ grpn_name1{b} '_' abs_band_name1{c}]);
         end
       end
       
       for b=1:length(grpn_name1);
            fprintf(fid,';%s',[ grpn_name1{b} '_abs_tot']);
       end
        
        for b=1:length(grpn_name1);
         for c=1:length(rel_band_name1);
            fprintf(fid,';%s',[ grpn_name1{b} '_' rel_band_name1{c}]);
         end
        end

        for b=1:length(grpn_name1);
            fprintf(fid,';%s',[ grpn_name1{b} '_SEF']);
        end
        
        fprintf(fid,';%s',[])
        
        fprintf(fid,'\n');
    
     fprintf(fid,';%s',inpfilename1);
     fprintf(fid,[repmat(';%f',[1 length(abs_pwr_mat1)])],abs_pwr_mat1);
     fprintf(fid,[repmat(';%f',[1 length(abs_power_total1)])],abs_pwrlog_total1);
     fprintf(fid,[repmat(';%f',[1 length(rel_pwr_mat1)])],rel_pwr_mat1);
     fprintf(fid,[repmat(';%f',[1 length(SEF1)])],SEF1);
     fprintf(fid,'\n');

      
end    
fclose(fid);

