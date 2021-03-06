%% Clears and Comments
clear;
%clc;
clf;

%Inputs: Donor lifetime.
%Outputs: Background, FRET lifetime, FRET and Non-FRET fractions
%V8: added w1min
%Added extract fitting dimension.

%Need to add:
%1) For window factor, don't use best in each dimension, but instead use best in parameter space

%2) Marginalize to find frac_est, then apply scale factor
tic
for jind = 1:1
    
    for expt = 1:1
        % fprintf('\njind is :%f expt is :%f\n', jind, expt);
        
        if expt == 1
            
          cyclesmax = 1;

            str(1).name = 'm0-spot1.sdt';
            irfname = 'irf-157pm.sdt';
            data_shift_name = 'm0-spot1.sdt';
        end
        

        
        
        for index = 1:cyclesmax
            
            %fprintf('\nindex is :%f    filename: is %s\n', index,str(index).name);
            %fprintf('\nexpt is :%f index is :%f\n', expt, index);
            
            
            %% Make IRF, wigs, IRF shift, wigs shift, and extract signal, read in data path
            
            sysinfo = 0; % Set to 1 if you want to force a rerun of make_irf_wig_ext
            pth_irf = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-6-13\'; %Type the file path of the Instrument Response Function (IRF) here.
            %irfname = 'irf-641pm.sdt';
            
            pth_wigs = 'Z:\bkaye\Bayes_2013\data\2014-2-13\';
            wigsname = 'wigs_1p5e6.sdt';
            
            pth_ext = 'Z:\bkaye\Bayes_2013\data\2013-5-16\';
            extname = 'm13_scan1.sdt';
            
            pth_data_for_shift = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-6-13\';
            %data_shift_name = 'm0-e1-s1-spot1.sdt';
            
            cppout = fopen('SysInfo.txt');
            [old_filenames, count] = fscanf(cppout, '%s');
            new_filenames = [pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name];
            fclose(cppout);
            
            if isequal(new_filenames,old_filenames) && sysinfo == 0;
                
            else
                sfrac = 0.2; shiftmin = -10; shiftmax = 2;
                w2step_shift = 0.01; w2min_shift = 3.5; w2max_shift = 4;
                backstep = 0.01; backmin = 0.03; backmax = 0.1;
                wigstep = 2; wigmin = 0; wigmax = 0;
                
                [~, ~, ~, ~, ~, ~, ~, ~, ~] = make_irf_wig_ext(pth_irf, irfname, pth_wigs, wigsname, pth_ext, extname, data_shift_name, pth_data_for_shift, sfrac, shiftmin, shiftmax, w2step_shift,w2min_shift, w2max_shift, backstep, backmin, backmax, wigstep, wigmin, wigmax);
                fileID = fopen('SysInfo.txt','w');
                fprintf(fileID,'%s%s\n%s%s\n%s%s\n%s%s\n',pth_irf, irfname, pth_wigs, wigsname, pth_ext,extname, pth_data_for_shift, data_shift_name);
                fclose(fileID);
                %fprintf('%s', data_shift_name);
            end
            
            %% Load in IRF, wigs
            
            tempf=load('Z:\bkaye\Bayes_2013\mat_files\current.mat','-mat','irf','bneed',...
                'pulsewb', 'irf', 'irfsim', 'irf_pdf', 'tmini', 'tmaxi', 'ext','irfname','wigsb','gab');
            irf=tempf(1).irf;        %Probability of delay (pdf_IRF)
            brem = tempf(1).bneed;   %Number of bins to remove from the 12.5ns to match BH data
            bins = tempf(1).pulsewb; %Number of bins that make up 12.5ns
            tmini = tempf(1).tmini;
            tmaxi = tempf(1).tmaxi;
            ext = tempf(1).ext;
            wig = tempf(1).wigsb;
            ga = tempf(1).gab;
            binskeep = bins-brem;
            

            %%
            
            if jind ==1
            pth_data = 'C:\Users\Bryan\Documents\MATLAB\Bayes_2014\data\2014-6-13\';
            dataname = str(index).name;
            thrt = 0; %0 means do not threshold, 1 means threshold , 2 means treshold for desired photon amount
            imthr = 3;%amount to threshol by. Thresh is max(int)/imthr
            DesPho = 200000;%jind*100000;
            
            [p] = read_imagethr_14ns(tmini,tmaxi,dataname,pth_data,imthr,thrt,DesPho);
            
            save('Z:\bkaye\cluster\data.mat', 'p');
%             datatemp = load('Z:\bkaye\cluster\data.mat','-mat','p');
%             p = datatemp(1).p;

            end
        end
    end
end
