function [filt_data_neu, filt_data_neu_divided,...
          filt_data_hep, filt_data_hep_divided] = revised_tech_interview_new(input_filename)
    %% Input Data
    % Import Filename 
    filename = input_filename;
    % Find Out How Many Spreadsheets
    [~, sheets] = xlsfinfo(filename);
        
    % Parameters Settings For Further Purposes
    % For Looping Purpose
    n_drugs = length(sheets);
    
    % For Plotting Purpose
    % setup figure counter
    pic_counter_Hep = 1;
    pic_counter_NEU = 3;
    
    %% 1.a Dose Response Curve (Including Both Original and Divided by Area Ones)
% %         ## this is not good, as i_drug=2 will erase i_drug=1 result
% %         ## more importantly, what if NEU have different number of Data IDs from HEP? 
% %         ## i.e., if Neurons have 3 datasets while HepG2 has 4 datasets?
% %             ## 3 dataset for Neuron, 3 dataset for HepG2, but you make zeros(2,6) : 12 spaces
% %             ## you might want to make zeros(2,3) because 2 for cell type, 3 for data ID
% %                 ## because we have 3 datasets for each cell type
% %         ## an example of solution: since did includes the information of cell type (not independent)
% %         ## mean_zero_array = zeros(n_dids,1) because each did will have only one value.
% %         ## if you want to save information of which did corresponds to which cell type, 
% %             ## you can use anothe rvariable cellType = strings(n_dids,1).
 
    % Data Reading Process
    for i_drug=1:n_drugs
        % Reading Data by Spreadsheet
        data_current = readtable(filename,'Sheet',i_drug);

        % Create Variable for Storing Predictors and Responses
        pred_hep = [];
        resp_hep = [];
        pred_neu = [];
        resp_neu = [];
        pred_hep_divided = [];
        resp_hep_divided = [];
        pred_neu_divided = [];
        resp_neu_divided = [];
        
        % Calculate the Number of Data IDs in This SpreadSheet
        % Number of Data IDs Calculation
        % Retrieve the data id for each sheet and calculate using numel
        % function
        data_in_ids = unique(data_current.did);
        n_dids = numel(data_in_ids);
        
        % Looping through Each Batch 
        for i_did=1:n_dids
           % Filter Out Unnecessary Data 
           % (preserve only the data in current batch)
           did_current = data_in_ids(i_did);
           data_current_did = data_current(string(data_current.did) == did_current,:);
           % Cell Type Indicator (For Later Storage Purpose) 
           % (1 for HepG2 and 2 for Neu in our case  but it can be changed by little modification in the previous portion)
           if contains(data_current_did.cell(~contains(data_current_did.cell, 'blank'),:), 'HepG2')
               pic_counter = pic_counter_Hep + 4*(i_drug-1);
               cellType_str = 'HepG2';
           elseif contains(data_current_did.cell(~contains(data_current_did.cell, 'blank'),:), 'neuro')
               pic_counter = pic_counter_NEU + 4*(i_drug-1);
               cellType_str = 'Neuro';
           end

           % Calculation Portion
           % First Step: Zero Concentration Data
           % Grab Out Only Zero Concentration Data
           % (For filtering out data that are blank or incubator cases)
           data_current_did_zero = data_current_did(data_current_did.conc==0 &...
                                                    ~contains(data_current_did.cell, 'blank') &...
                                                    ~contains(data_current_did.cell, '-inc'),:);
           % Original Data Portion
           data_current_did_zero_lumin = data_current_did_zero.lumin;
           % Calculate the Mean and CI_95 for Zero Concentration Data
           % (Make sure it's divided by itself to make it percentile)
           ci = 0.95 ;
           a = 1 - ci;
           n_zero = numel(data_current_did_zero_lumin);
           Tz_multiplier = tinv(1-a/2, n_zero-1);
           ci_95_current_did_zero = Tz_multiplier*std(data_current_did_zero_lumin)/sqrt(n_zero);
           mean_current_did_zero_temp = mean(data_current_did_zero_lumin);
           ci_95_current_did_zero = ci_95_current_did_zero/mean_current_did_zero_temp;
           mean_current_did_zero = mean_current_did_zero_temp/mean_current_did_zero_temp;
           % Plot Current Zero Concentration 
           figure(pic_counter);
           title(['Dose vs Concentration Plot for ', cellType_str, ' (Original Version)' ]);
           clr = lines(length(data_in_ids));
           % Plot the data
           % Set up a fake x axis for zero concentration by finding the
           % minimum value of all nonzero concentration in the spread sheet and divided
           % it by 10^6
           hold on
           factor = 1e6;
           conc_min = min(data_current.conc(data_current.conc>0,:));
           fake_x_zero = conc_min/factor;
           e = errorbar(fake_x_zero,mean_current_did_zero,ci_95_current_did_zero,'color',clr(i_did,:), 'linestyle','none');
           e.LineStyle = 'none';
           % Plot display setting
           set(gca, 'XScale', 'log')
           set(gca, 'YScale', 'log')
           % Label the plot
           xlabel('Concentration')
           ylabel('Luminance (%)')
           grid on
           % Store the Data into Predictor and Response Vector into
           % Corresponding Vector
           if strcmp(cellType_str,'HepG2')
               if isempty(pred_hep) && isempty(resp_hep)
                  pred_hep = fake_x_zero;
                  resp_hep = mean_current_did_zero;
               elseif ~isempty(pred_hep) && ~isempty(resp_hep)
                  pred_hep = cat(1,pred_hep,fake_x_zero);
                  resp_hep = cat(1,resp_hep,mean_current_did_zero);
               else
                   error('Something goes wrong in the HepG2 Zero Concentration Matrices')
               end
           elseif strcmp(cellType_str,'Neuro')
               if isempty(pred_neu) && isempty(resp_neu)
                  pred_neu = fake_x_zero;
                  resp_neu = mean_current_did_zero;
               elseif ~isempty(pred_neu) && ~isempty(resp_neu)
                  pred_neu = cat(1,pred_neu,fake_x_zero);
                  resp_neu = cat(1,resp_neu,mean_current_did_zero);
               else
                   error('Something goes wrong in the Neuro Zero Concentration Matrices')
               end
           else
               error('The function needs to be modified. The celltype is changed.')
           end
           % Divided by Area Portion
           % Grab out Area Data
           data_current_did_zero_area = data_current_did_zero.area;
           data_current_did_zero_lumin_divided = data_current_did_zero_lumin./data_current_did_zero_area;
           ci_95_current_did_zero_divided = Tz_multiplier*std(data_current_did_zero_lumin_divided)/sqrt(n_zero);
           mean_current_did_zero_divided_temp = mean(data_current_did_zero_lumin_divided);
           ci_95_current_did_zero_divided = ci_95_current_did_zero_divided/mean_current_did_zero_divided_temp;
           mean_current_did_zero_divided = mean_current_did_zero_divided_temp/mean_current_did_zero_divided_temp;
           % Plot Current Zero Concentration Divided by Area
           figure(pic_counter + 1);
           title(['Dose vs Concentration Plot for ', cellType_str, ' (Divided by Area)' ]);
           hold on
           e = errorbar(fake_x_zero,mean_current_did_zero_divided,ci_95_current_did_zero_divided,'color',clr(i_did,:), 'linestyle','none');
           e.LineStyle = 'none';
           % Plot display setting
           set(gca, 'XScale', 'log')
           set(gca, 'YScale', 'log')
           % Label the plot
           xlabel('Concentration')
           ylabel('Luminance (%)')
           grid on
           % Store the Data into Predictor and Response Vector into
           % Corresponding vector
           if strcmp(cellType_str,'HepG2')
               if isempty(pred_hep_divided) && isempty(resp_hep_divided)
                   pred_hep_divided = fake_x_zero;
                   resp_hep_divided = mean_current_did_zero_divided;
               elseif ~isempty(pred_hep_divided) && ~isempty(resp_hep_divided)
                   pred_hep_divided = cat(1,pred_hep_divided,fake_x_zero);
                   resp_hep_divided = cat(1,resp_hep_divided,mean_current_did_zero_divided);
               else
                   error('Something goes wrong in the HepG2 Zero Concentration Divided Matrices')
               end
           elseif strcmp(cellType_str,'Neuro')
               if isempty(pred_neu_divided) && isempty(resp_neu_divided)
                   pred_neu_divided = fake_x_zero;
                   resp_neu_divided = mean_current_did_zero_divided;
               elseif ~isempty(pred_neu_divided) && ~isempty(resp_neu_divided)
                   pred_neu_divided = cat(1,pred_neu_divided,fake_x_zero);
                   resp_neu_divided = cat(1,resp_neu_divided,mean_current_did_zero_divided);   
               else
                   error('Something goes wrong in the Neuro Zero Concentration Divided Matrices')
               end
           else
               error('The function needs to be modified. The celltype is changed.')
           end
           
           
           % Second Step: Incubator Data
           % Grab Out Only Incubator Data
           % (For filtering out data that are zero concentration or blank cases)
           data_current_did_inc = data_current_did(contains(data_current_did.cell, '-inc'),:);
           % Original Data Portion
           data_current_did_inc_lumin = data_current_did_inc.lumin;
           % Calculate the Mean and CI_95 for Blank Data
           % (Make sure it's divided by zero concentration data to make it percentile)
           mean_current_did_inc = mean(data_current_did_inc_lumin);
           mean_current_did_inc = mean_current_did_inc/mean_current_did_zero_temp;
           ci = 0.95 ;
           a = 1 - ci;
           n_inc = numel(data_current_did_inc_lumin);
           Tinc_multiplier = tinv(1-a/2, n_inc-1);
           ci_95_current_did_inc = Tinc_multiplier*std(data_current_did_inc_lumin)/sqrt(n_inc);
           ci_95_current_did_inc = ci_95_current_did_inc/mean_current_did_zero_temp;
           % Plotting Current Incubator Data
           figure(pic_counter);
           clr = lines(length(data_in_ids));
           % Set up a fake x axis for zero concentration by finding the
           % minimum value of all nonzero concentration in the spread sheet and divided
           % it by 10^7
           hold on
           factor = 1e7;
           conc_min = min(data_current.conc(data_current.conc>0,:));
           fake_x_inc = conc_min/factor;
           e = errorbar(fake_x_inc,mean_current_did_inc,ci_95_current_did_inc,'color',clr(i_did,:), 'linestyle','none');
           e.LineStyle = 'none';
           % Plot display setting
           set(gca, 'XScale', 'log')
           set(gca, 'YScale', 'log')
           % Label the plot
           xlabel('Concentration')
           ylabel('Luminance (%)')
           grid on
           % Store the Data into Predictor and Response Vector into
           % Corresponding Vecotr
           if strcmp(cellType_str,'HepG2')
              pred_hep = cat(1,pred_hep,fake_x_inc);
              resp_hep = cat(1,resp_hep,mean_current_did_inc);
           elseif strcmp(cellType_str,'Neuro')
               pred_neu = cat(1,pred_neu,fake_x_inc);
               resp_neu = cat(1,resp_neu,mean_current_did_inc);
           else
               error('The function needs to be modified. The celltype is changed.')
           end
               
           % Divided by Area Portion
           % Grab out Area Data
           data_current_did_inc_area = data_current_did_inc.area;
           data_current_did_inc_lumin_divided = data_current_did_inc_lumin./data_current_did_inc_area;
           mean_current_did_inc_divided = mean(data_current_did_inc_lumin_divided);
           mean_current_did_inc_divided = mean_current_did_inc_divided/mean_current_did_zero_divided_temp;
           ci_95_current_did_inc_divided = Tinc_multiplier*std(data_current_did_inc_lumin_divided)/sqrt(n_inc);
           ci_95_current_did_inc_divided = ci_95_current_did_inc_divided/mean_current_did_zero_divided_temp;
           % Plot Current Incubator Data Divided by Area 
           figure(pic_counter+1);
           e = errorbar(fake_x_inc,mean_current_did_inc_divided,ci_95_current_did_inc_divided,'color',clr(i_did,:), 'linestyle','none');
           e.LineStyle = 'none';
           % Plot display setting
           set(gca, 'XScale', 'log')
           set(gca, 'YScale', 'log')
           % Label the plot
           xlabel('Concentration')
           ylabel('Luminance (%)')
           grid on
           % Store the Data into Predictor and Response Vector into
           % Corresponding Vector
           if strcmp(cellType_str,'HepG2')
              pred_hep_divided = cat(1,pred_hep_divided,fake_x_inc);
              resp_hep_divided = cat(1,resp_hep_divided,mean_current_did_inc_divided);
           elseif strcmp(cellType_str,'Neuro')
              pred_neu_divided = cat(1,pred_neu_divided,fake_x_inc);
              resp_neu_divided = cat(1,resp_neu_divided,mean_current_did_inc_divided);
           else
               error('The function need to be modified. The celltype is changed.')
           end

           % Third Step: Concentration Data
           % Grab Out Only Concentration Data 
           % (With Concentrationos That Are Nonzero)
           data_current_did_CONCS = data_current_did(data_current_did.conc>0,:);
           % Detect How Many Concentrations Were Tested in This Batch
           curr_concs = unique(data_current_did_CONCS.conc); 
           n_concs = numel(curr_concs);
           % Looping through each concentration
           for i_conc=1:n_concs
              % Grab Out the Current Concentration Data
              data_current_did_currconc = data_current_did_CONCS(data_current_did_CONCS.conc==curr_concs(i_conc),:);
              % Original Data Portion
              data_current_did_currconc_lumin = data_current_did_currconc.lumin;
              % Calculate the Mean and CI_95 for Current Concentration Data
              % (Make sure it's divided by zero concentration data to make it percentile)
              mean_current_did_currconc = mean(data_current_did_currconc_lumin);
              mean_current_did_currconc = mean_current_did_currconc/mean_current_did_zero_temp;
              ci = 0.95 ;
              a = 1 - ci;
              n_currconc = numel(data_current_did_currconc_lumin);
              Tcurrconc_multiplier = tinv(1-a/2, n_currconc-1);
              ci_95_current_did_currconc = Tcurrconc_multiplier*std(data_current_did_currconc_lumin)/sqrt(n_currconc);
              ci_95_current_did_currconc = ci_95_current_did_currconc/mean_current_did_zero_temp;
              % Plot Current Concentration Data
              current_conc = curr_concs(i_conc);
              clr = lines(length(data_in_ids));
              figure(pic_counter);
              hold on
              e = errorbar(current_conc,mean_current_did_currconc,ci_95_current_did_currconc,'color',clr(i_did,:), 'linestyle','none');
              e.LineStyle = 'none';
              % plot display setting
              set(gca, 'XScale', 'log')
              set(gca, 'YScale', 'log')
              % label the plot
              xlabel('Concentration')
              ylabel('Luminance (%)')
              grid on
              % Store the Data into Predictor and Response Vector into
              % Corresponding Vector
              if strcmp(cellType_str,'HepG2')
                 pred_hep = cat(1,pred_hep,current_conc);
                 resp_hep = cat(1,resp_hep,mean_current_did_currconc);
              elseif strcmp(cellType_str,'Neuro')
                 pred_neu = cat(1,pred_neu,current_conc);
                 resp_neu = cat(1,resp_neu,mean_current_did_currconc);
              else
                  error('The function need to be modified. The celltype is changed.')
              end
              
              % Divided by Area Portion
              data_current_did_currconc_area = data_current_did_currconc.area;
              data_current_did_currconc_lumin_divided = data_current_did_currconc_lumin./data_current_did_currconc_area;
              mean_current_did_currconc_divided = mean(data_current_did_currconc_lumin_divided);
              mean_current_did_currconc_divided = mean_current_did_currconc_divided/mean_current_did_zero_divided_temp;
              ci_95_current_did_currconc_divided = Tcurrconc_multiplier*std(data_current_did_currconc_lumin_divided)/sqrt(n_currconc);
              ci_95_current_did_currconc_divided = ci_95_current_did_currconc_divided/mean_current_did_zero_divided_temp;
              % Plot Divided Data
              figure(pic_counter+1);
              hold on
              e = errorbar(current_conc,mean_current_did_currconc_divided,ci_95_current_did_currconc_divided,'color',clr(i_did,:), 'linestyle','none');
              e.LineStyle = 'none';
              % plot display setting
              set(gca, 'XScale', 'log')
              set(gca, 'YScale', 'log')
              % label the plot
              xlabel('Concentration')
              ylabel('Luminance (%)')
              grid on
              % Store the Data into Predictor and Response Vector into
              % Corresponding Vector 
              if strcmp(cellType_str,'HepG2')
                 pred_hep_divided = cat(1,pred_hep_divided,current_conc);
                 resp_hep_divided = cat(1,resp_hep_divided,mean_current_did_currconc_divided);
              elseif strcmp(cellType_str,'Neuro')
                 pred_neu_divided = cat(1,pred_neu_divided,current_conc);
                 resp_neu_divided = cat(1,resp_neu_divided,mean_current_did_currconc_divided);
              else
                  error('The function need to be modified. The celltype is changed.')
              end
           end      
        end
        
        % Forth Step: Blank Data
        % HepG2 Case
        % Calculate the Mean and CI95 of HepG2 CellType 
        % (Make sure it's divided by the overall zero concentration mean)
        data_current_Hep_zero = data_current(string(data_current.cell) == 'HepG2'&...
                                             string(data_current.cell) ~= 'blank'&...
                                            ~contains(data_current.cell,'-inc')&...
                                             data_current.conc==0,:);
        data_current_Hep_zero_lumin = data_current_Hep_zero.lumin;
        mean_current_Hep_zero = mean(data_current_Hep_zero_lumin);
        % Grab Out All HepG2 Batches Blank Data
        data_current_Hep_bla = data_current(contains(data_current.did, 'HEP')&...
                                            string(data_current.cell)=='blank',:);
        % Original Data Portion
        data_current_Hep_bla_lumin = data_current_Hep_bla.lumin;
        % Setup Celltype Indicator (For HepG2, 1)
        mean_current_Hep_bla = mean(data_current_Hep_bla_lumin);
        mean_current_Hep_bla = mean_current_Hep_bla/mean_current_Hep_zero;
        ci = 0.95 ;
        a = 1 - ci;
        n_bla = numel(data_current_Hep_bla_lumin);
        Tbla_multiplier = tinv(1-a/2, n_bla-1);
        ci_95_current_Hep_bla = Tbla_multiplier*std(data_current_Hep_bla_lumin)/sqrt(n_bla);
        ci_95_current_Hep_bla = ci_95_current_Hep_bla/mean_current_Hep_zero;
        % Plot HepG2 Blank Data
        figure(pic_counter_Hep + 4*(i_drug-1));
        hold on
        upp = (mean_current_Hep_bla + ci_95_current_Hep_bla);
        yline(upp, '--')
        hold on
        dow = (mean_current_Hep_bla - ci_95_current_Hep_bla);
        yline(dow, '--')
        hold on
        

        % Model Fitting 
        % Construct Models for Fitting
        % Setting Up Fixed-effect Coefficients
        % Scenario 1: x0 only for random effect

        % Preprocessing the Data
        % Get All Unique Concentration Data
        unique_pred_hep = unique(pred_hep);
        % Create New Storage Variable for Storing Filtered Data
        filt_unique_pred_hep = [];
        filt_unique_resp_hep = [];
        % Looping through Each Unique Concentration
        num_uniq_pred_hep = size(unique_pred_hep,1);
        for i_uniq_pred_hep=1:num_uniq_pred_hep
            curr_unique_pred_hep = unique_pred_hep(i_uniq_pred_hep);
            % Concatenate the current concentration into the new storage
            % variable for predictor
            filt_unique_pred_hep = cat(1,filt_unique_pred_hep,curr_unique_pred_hep);
            % Grab out all indices that are corresponding to this
            % concentration 
            curr_unique_pred_hep_ind = find(pred_hep==curr_unique_pred_hep);
            % Get the frequency of the corresponding concentration
            num_curr_unique_pred_hep = size(curr_unique_pred_hep_ind,1);
            % If there's more than one, average through all the
            % corresponding data and store it into the storage variable
            if num_curr_unique_pred_hep > 1
               avg_curr_unique_hep_resp = mean(resp_hep(curr_unique_pred_hep_ind));
               filt_unique_resp_hep = cat(1,filt_unique_resp_hep, avg_curr_unique_hep_resp);
            % If the corresponding concentration only have one
            % corresponding response, just store it directly.
            else
               filt_unique_resp_hep = cat(1,filt_unique_resp_hep, resp_hep(curr_unique_pred_hep_ind));
            end
        end

        % Looping through each of the elements to check the longest
        % decreasing interval
        % Concatenate Predictor and Response Array Together to Form
        % Corresponding Predictor-Response Pairs Data
        filt_data_hep = cat(2, filt_unique_pred_hep, filt_unique_resp_hep);
        % Sort Out Data Array in Ascending Order of the First Column (Predictor Array)
        % and Grab Out Its Size
        sorted_filt_data_hep = sortrows(filt_data_hep);
        n_col_filt_data_hep = size(filt_data_hep,1);
        n_row_filt_data_hep = size(filt_data_hep,2);
        % Create Variable for Storing Maximum Interval of Response
        % Descending and Store the Corresponding Start and End predictors 
        max_interval = 0;
        max_start_point = 0;
        max_end_point = 0;
        % Create Variable to Grab Out the Starting and Ending Point of the Current Interval
        Start_point = 0;
        End_point = 0;
        % Create Variable for Storing Starting and Ending Index (for
        % pointing purpose at the data matrix)
        Start_ind = 0;
        End_ind = 0;
        % Loop through the data
        for i_col_filt_data_hep = 1:n_col_filt_data_hep
            % Check if the interval is descending
            % If this is the first looping, save the corresponding
            % predictor value and index and skip the investigation
            if i_col_filt_data_hep == 1
                Start_point = sorted_filt_data_hep(i_col_filt_data_hep,1);
                Start_ind = i_col_filt_data_hep;
            % If this is not the first looping, check if the interval is
            % descending
            else
                % Grab out current tested response and the starting point
                curr_resp = sorted_filt_data_hep(i_col_filt_data_hep,n_row_filt_data_hep);
                start_resp = sorted_filt_data_hep(Start_ind,n_row_filt_data_hep);
                % If the current interval is descending, update the current
                % point to endpoint and update current index to end index
                if curr_resp < start_resp && i_col_filt_data_hep ~= n_col_filt_data_hep
                    End_point = sorted_filt_data_hep(i_col_filt_data_hep,1);
                    End_ind = i_col_filt_data_hep;
                % If the current interval is not descending (or anymore), 
                % compare the current interval with existed maximum
                % interval.
                else
                    % Calculate the current interval
                    curr_interval = End_point-Start_point;
                    % Compare it with the current existed maximum interval
                    % If it's larger, update the existed maximum interval
                    % value.
                    if curr_interval > max_interval || (curr_resp < start_resp && i_col_filt_data_hep == n_col_filt_data_hep)
                        max_interval = curr_interval;
                        max_start_point = Start_point;
                        max_end_point = End_point;
                        if i_col_filt_data_hep == n_col_filt_data_hep
                            guess_x0 = max_end_point;
                        else
                            guess_x0 = (max_start_point+max_end_point)/2;
                        end
                    end
                    Start_point = sorted_filt_data_hep(i_col_filt_data_hep,1);
                    Start_ind = i_col_filt_data_hep;
                    End_point = sorted_filt_data_hep(i_col_filt_data_hep,1);
                    End_ind = i_col_filt_data_hep;
                end
            end
        end
        lumin_high = max(filt_unique_resp_hep);
        lumin_low = min(filt_unique_resp_hep);
        % Sample a Point in Data (But Avoid Choosing lumin_high
        % and lumin_low)
        lumin_for_sample = filt_unique_resp_hep(filt_unique_resp_hep~=lumin_high & filt_unique_resp_hep~=lumin_low);
        lumin_for_sample = sort(lumin_for_sample);
        length_lumin_sample = size(lumin_for_sample,1);
        lumin_samp = lumin_for_sample(ceil(length_lumin_sample/2));
        [sam_ind_x,sam_ind_y] = find(filt_unique_resp_hep==lumin_samp);
        conc_sam = filt_unique_pred_hep(sam_ind_x,sam_ind_y);
        if mean_current_Hep_bla > lumin_low
            cal_var1 = log((lumin_high-mean_current_Hep_bla)/(lumin_samp-mean_current_Hep_bla)-1);
        else
            cal_var1 = log((lumin_high-lumin_low)/(lumin_samp-lumin_low)-1);
        end
        cal_var2 = (-1)/(conc_sam-guess_x0);
        k = cal_var1*cal_var2;
        % Determine the minimum (lower) boundary of the function
        if mean_current_Hep_bla > lumin_low
            % Upper Bound Setting
            L = lumin_high-mean_current_Hep_bla; % since the upper bound is supposed to be 100%  
            model = @(phi,c)(L./(1+exp(-k.*(c-phi)))+mean_current_Hep_bla);
        else
            % Upper Bound Setting
            L = lumin_high-lumin_low; % since the upper bound is supposed to be 100%  
            model = @(phi,c)(L./(1+exp(-k.*(c-phi)))+lumin_low);
        end
        % initial guess
        phi0 = guess_x0; 
        [phi_Hep,res_Hep] = nlinfit(pred_hep,resp_hep,model,phi0);
        % Plotting Outcomes
        min_pre = min(pred_hep);
        max_pre = max(pred_hep);
        sim_input = min_pre:0.00001:max_pre;
        % Plotting the model
        figure(pic_counter_Hep + 4*(i_drug-1));
        hold on
        plot(sim_input,model(phi_Hep,sim_input),'k','LineWidth',2)
        out_hep = sim_input;model(phi_Hep,sim_input);
        % plot display setting
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        hold on
        % Outputing Mean Squared Error
        num0bs = length(pred_hep);
        numParams = 1;
        df = num0bs - numParams;
        mse = (res_Hep'*res_Hep)/df;
        disp('The mean squared error of the first model on HepG2 is:')
        disp(mse)
        
        % Divided by Area Portion
        % In addition to the area of current celltype, the mean of divided zero
        % concentration data is also needed
        % To get that, we first need to grab out all zero concentration in
        % this celltype and then divided its luminscences by its areas 
        data_current_Hep_zero_lumin = data_current_Hep_zero.lumin;
        data_current_Hep_zero_area = data_current_Hep_zero.area;
        data_current_Hep_zero_divided = data_current_Hep_zero_lumin./data_current_Hep_zero_area;
        % After that, we calculate its mean
        mean_current_Hep_zero_divided = mean(data_current_Hep_zero_divided);
        % For blank cases, it needs to be divided by the mean area value of
        % the current celltype.
        % We then first grab out the area data and calculate its mean
        data_current_Hep_area = data_current.area(contains(data_current.cell,'HepG2'),:);
        mean_area_current_Hep = mean(nonzeros(data_current_Hep_area));
        % After doing so, we then divided the blank data with the mean area
        % of this cellType
        data_current_Hep_bla_lumin_divided = data_current_Hep_bla_lumin/mean_area_current_Hep;
        mean_current_Hep_bla_divided = mean(data_current_Hep_bla_lumin_divided);
        % To achieve the final percentile data, we then utilize the mean of
        % overall divided zero concentration and divided our obtained mean blank value by it 
        mean_current_Hep_bla_divided = mean_current_Hep_bla_divided/mean_current_Hep_zero_divided;
        ci_95_current_Hep_bla_divided = Tbla_multiplier*std(data_current_Hep_bla_lumin_divided)/sqrt(n_bla);
        ci_95_current_Hep_bla_divided = ci_95_current_Hep_bla_divided/mean_current_Hep_zero_divided;
        % Plot HepG2 Data Divided by Area
        figure(pic_counter_Hep + 4*(i_drug-1)+1);
        hold on
        upp = (mean_current_Hep_bla_divided + ci_95_current_Hep_bla_divided);
        yline(upp, '--')
        hold on
        dow = (mean_current_Hep_bla_divided - ci_95_current_Hep_bla_divided);
        yline(dow, '--')
        hold on
        

        % Model Fitting 
        % Construct Models for Fitting
        % Setting Up Fixed-effect Coefficients 
        % Scenario 1: x0 only for random effect

        % Preprocessing the Data
        % Get All Unique Concentration Data
        unique_pred_hep_divided = unique(pred_hep_divided);
        % Create New Storage Variable for Storing Filtered Data
        filt_unique_pred_hep_divided = [];
        filt_unique_resp_hep_divided = [];
        % Looping through Each Unique Concentration
        num_uniq_pred_hep_divided = size(unique_pred_hep_divided,1);
        for i_uniq_pred_hep_divided=1:num_uniq_pred_hep_divided
            curr_unique_pred_hep_divided = unique_pred_hep_divided(i_uniq_pred_hep_divided);
            % Concatenate the current concentration into the new storage
            % variable for predictor
            filt_unique_pred_hep_divided = cat(1,filt_unique_pred_hep_divided,curr_unique_pred_hep_divided);
            % Grab out all indices that are corresponding to this
            % concentration 
            curr_unique_pred_hep_divided_ind = find(pred_hep_divided==curr_unique_pred_hep_divided);
            % Get the frequency of the corresponding concentration
            num_curr_unique_pred_hep_divided = size(curr_unique_pred_hep_divided_ind,1);
            % If there's more than one, average through all the
            % corresponding data and store it into the storage variable
            if num_curr_unique_pred_hep_divided > 1
               avg_curr_unique_hep_resp_divided = mean(resp_hep_divided(curr_unique_pred_hep_divided_ind));
               filt_unique_resp_hep_divided = cat(1,filt_unique_resp_hep_divided, avg_curr_unique_hep_resp_divided);
            % If the corresponding concentration only have one
            % corresponding response, just store it directly.
            else
               filt_unique_resp_hep_divided = cat(1,filt_unique_resp_hep_divided, resp_hep_divided(curr_unique_pred_hep_divided_ind));
            end
        end

        % Looping through each of the elements to check the longest
        % decreasing interval
        % Concatenate Predictor and Response Array Together to Form
        % Corresponding Predictor-Response Pairs Data
        filt_data_hep_divided = cat(2, filt_unique_pred_hep_divided, filt_unique_resp_hep_divided);
        % Sort Out Data Array in Ascending Order of the First Column (Predictor Array)
        % and Grab Out Its Size
        sorted_filt_data_hep_divided = sortrows(filt_data_hep_divided);
        n_col_filt_data_hep_divided = size(sorted_filt_data_hep_divided,1);
        n_row_filt_data_hep_divided = size(sorted_filt_data_hep_divided,2);
        % Create Variable for Storing Maximum Interval of Response
        % Descending and Store the Corresponding Start and End predictors 
        max_interval = 0;
        max_start_point = 0;
        max_end_point = 0;
        % Create Variable to Grab Out the Starting and Ending Point of the Current Interval
        Start_point = 0;
        End_point = 0;
        % Create Variable for Storing Starting and Ending Index (for
        % pointing purpose at the data matrix)
        Start_ind = 0;
        End_ind = 0;
        % Loop through the data
        for i_col_filt_data_hep_divided = 1:n_col_filt_data_hep_divided
            % Check if the interval is descending
            % If this is the first looping, save the corresponding
            % predictor value and index and skip the investigation
            if i_col_filt_data_hep_divided == 1
                Start_point = sorted_filt_data_hep_divided(i_col_filt_data_hep_divided,1);
                Start_ind = i_col_filt_data_hep_divided;
            % If this is not the first looping, check if the interval is
            % descending
            else
                % Grab out current tested response and the starting point
                curr_resp = sorted_filt_data_hep_divided(i_col_filt_data_hep_divided,n_row_filt_data_hep_divided);
                start_resp = sorted_filt_data_hep_divided(Start_ind,n_row_filt_data_hep_divided);
                % If the current interval is descending, update the current
                % point to endpoint and update current index to end index
                if curr_resp < start_resp && i_col_filt_data_hep_divided~=n_col_filt_data_hep_divided
                    End_point = sorted_filt_data_hep_divided(i_col_filt_data_hep_divided,1);
                    End_ind = i_col_filt_data_hep_divided;
                % If the current interval is not descending (or anymore), 
                % compare the current interval with existed maximum
                % interval.
                else
                    % Calculate the current interval
                    curr_interval = End_point-Start_point;
                    % Compare it with the current existed maximum interval
                    % If it's larger, update the existed maximum interval
                    % value.
                    if curr_interval > max_interval || (curr_resp < start_resp && i_col_filt_data_hep_divided==n_col_filt_data_hep_divided)
                        max_interval = curr_interval;
                        max_start_point = Start_point;
                        max_end_point = End_point;
                        if i_col_filt_data_hep_divided==n_col_filt_data_hep_divided
                            guess_x0 = max_end_point;
                        else
                            guess_x0 = (max_start_point+max_end_point)/2;
                        end
                    end
                    Start_point = sorted_filt_data_hep_divided(i_col_filt_data_hep_divided,1);
                    Start_ind = i_col_filt_data_hep_divided;
                    End_point = sorted_filt_data_hep_divided(i_col_filt_data_hep_divided,1);
                    End_ind = i_col_filt_data_hep_divided;
                end
            end
        end
        lumin_high = max(filt_unique_resp_hep_divided);
        lumin_low = min(filt_unique_resp_hep_divided);
        % Randomly Sample a Point in Data (But Avoid Choosing lumin_high
        % and lumin_low)
        lumin_for_sample = filt_unique_resp_hep_divided(filt_unique_resp_hep_divided~=lumin_high & filt_unique_resp_hep_divided~=lumin_low);
        lumin_for_sample = sort(lumin_for_sample);
        length_lumin_sample = size(lumin_for_sample,1);
        lumin_samp = lumin_for_sample(ceil(length_lumin_sample/2));
        [sam_ind_x,sam_ind_y] = find(filt_unique_resp_hep_divided==lumin_samp);
        conc_sam = filt_unique_pred_hep_divided(sam_ind_x,sam_ind_y);
        if mean_current_Hep_bla_divided > lumin_low
            cal_var1 = log((lumin_high-mean_current_Hep_bla_divided)/(lumin_samp-mean_current_Hep_bla_divided)-1);
        else
            cal_var1 = log((lumin_high-lumin_low)/(lumin_samp-lumin_low)-1);
        end
        cal_var2 = (-1)/(conc_sam-guess_x0);
        k = cal_var1*cal_var2;
        % Determine the minimum (lower) boundary of the function
        if mean_current_Hep_bla_divided > lumin_low
            % Upper Bound Setting
            L = lumin_high-mean_current_Hep_bla_divided; % since the upper bound is supposed to be 100% 
            model2 = @(phi,c)(L./(1+exp(-k.*(c-phi)))+mean_current_Hep_bla_divided);
        else
            % Upper Bound Setting
            L = lumin_high-lumin_low; % since the upper bound is supposed to be 100% 
            model2 = @(phi,c)(L./(1+exp(-k.*(c-phi)))+lumin_low);
        end
        % initial guess
        phi0 = guess_x0; 
        [phi_Hep_divided,res_Hep_divided] = nlinfit(pred_hep_divided,resp_hep_divided,model2,phi0);
        % Plotting Outcomes
        min_pre = min(pred_hep_divided);
        max_pre = max(pred_hep_divided);
        sim_input = min_pre:0.00001:max_pre;
        % Plotting the model
        figure(pic_counter_Hep + 4*(i_drug-1)+1);
        hold on
        plot(sim_input,model2(phi_Hep_divided,sim_input),'k','LineWidth',2)
        out_hep_divided = model2(phi_Hep_divided,sim_input);
        % plot display setting
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        hold on
        % Outputing Mean Squared Error
        num0bs = length(pred_hep_divided);
        numParams = 1;
        df = num0bs - numParams;
        mse = (res_Hep_divided'*res_Hep_divided)/df;
        disp('The mean squared error of the first model on HepG2 (Divided) is:')
        disp(mse)
        
        % Neuro Cases
        data_current_NEU_bla_overall = data_current(contains(data_current.did, 'NEU')&...
                                                    contains(data_current.cell, 'blank'),:);
        data_current_NEU_bla_lumin = data_current_NEU_bla_overall.lumin;
        % Calculate the Mean and CI95 of HepG2 CellType 
        % (Make sure it's divided by the overall zero concentration mean)
        data_current_NEU_zero = data_current(contains(data_current.cell,'neuro')&...
                                             ~contains(data_current.cell, 'blank')&...
                                             data_current.conc==0,:);
        % Original Data Portion
        data_current_NEU_zero_lumin = data_current_NEU_zero.lumin;
        mean_current_NEU_zero = mean(data_current_NEU_zero_lumin);
        mean_current_NEU_bla = mean(data_current_NEU_bla_lumin);
        mean_current_NEU_bla = mean_current_NEU_bla/mean_current_NEU_zero;
        ci = 0.95 ;
        a = 1 - ci;
        n_bla = numel(data_current_NEU_bla_lumin);
        Tbla_multiplier = tinv(1-a/2, n_bla-1);
        ci_95_current_NEU_bla = Tbla_multiplier*std(data_current_NEU_bla_lumin)/sqrt(n_bla);
        ci_95_current_NEU_bla = ci_95_current_NEU_bla/mean_current_NEU_zero;
        % Plotting Neuro Plot
        figure(pic_counter_NEU + 4*(i_drug-1));
        hold on
        upp = (mean_current_NEU_bla + ci_95_current_NEU_bla);
        yline(upp, '--')
        hold on
        dow = (mean_current_NEU_bla - ci_95_current_NEU_bla);
        yline(dow, '--')
        hold on
        

        % Model Fitting 
        % Construct Models for Fitting
        % Setting Up Fixed-effect Coefficients
        % Scenario 1: x0 only for random effect

        % Preprocessing the Data
        % Get All Unique Concentration Data
        unique_pred_neu = unique(pred_neu);
        % Create New Storage Variable for Storing Filtered Data
        filt_unique_pred_neu = [];
        filt_unique_resp_neu = [];
        % Looping through Each Unique Concentration
        num_uniq_pred_neu = size(unique_pred_neu,1);
        for i_uniq_pred_neu=1:num_uniq_pred_neu
            curr_unique_pred_neu = unique_pred_neu(i_uniq_pred_neu);
            % Concatenate the current concentration into the new storage
            % variable for predictor
            filt_unique_pred_neu = cat(1,filt_unique_pred_neu,curr_unique_pred_neu);
            % Grab out all indices that are corresponding to this
            % concentration 
            curr_unique_pred_neu_ind = find(pred_neu==curr_unique_pred_neu);
            % Get the frequency of the corresponding concentration
            num_curr_unique_pred_neu = size(curr_unique_pred_neu_ind,1);
            % If there's more than one, average through all the
            % corresponding data and store it into the storage variable
            if num_curr_unique_pred_neu > 1
               avg_curr_unique_neu_resp = mean(resp_neu(curr_unique_pred_neu_ind));
               filt_unique_resp_neu = cat(1,filt_unique_resp_neu, avg_curr_unique_neu_resp);
            % If the corresponding concentration only have one
            % corresponding response, just store it directly.
            else
               filt_unique_resp_neu = cat(1,filt_unique_resp_neu, resp_neu(curr_unique_pred_neu_ind));
            end
        end

        % Looping through each of the elements to check the longest
        % decreasing interval
        % Concatenate Predictor and Response Array Together to Form
        % Corresponding Predictor-Response Pairs Data
        filt_data_neu = cat(2, filt_unique_pred_neu, filt_unique_resp_neu);
        % Sort Out Data Array in Ascending Order of the First Column (Predictor Array)
        % and Grab Out Its Size
        sorted_filt_data_neu = sortrows(filt_data_neu);
        n_col_sorted_filt_data_neu = size(sorted_filt_data_neu,1);
        n_row_sorted_filt_data_neu = size(sorted_filt_data_neu,2);
        % Create Variable for Storing Maximum Interval of Response
        % Descending and Store the Corresponding Start and End predictors 
        max_interval = 0;
        max_start_point = 0;
        max_end_point = 0;
        % Create Variable to Grab Out the Starting and Ending Point of the Current Interval
        Start_point = 0;
        End_point = 0;
        % Create Variable for Storing Starting and Ending Index (for
        % pointing purpose at the data matrix)
        Start_ind = 0;
        End_ind = 0;
        % Loop through the data
        for i_col_sorted_filt_data_neu = 1:n_col_sorted_filt_data_neu
            % Check if the interval is descending
            % If this is the first looping, save the corresponding
            % predictor value and index and skip the investigation
            if i_col_sorted_filt_data_neu == 1 
                Start_point = sorted_filt_data_neu(i_col_sorted_filt_data_neu,1);
                Start_ind = i_col_sorted_filt_data_neu;
            % If this is not the first looping, check if the interval is
            % descending
            else
                % Grab out current tested response and the starting point
                curr_resp = sorted_filt_data_neu(i_col_sorted_filt_data_neu,n_row_sorted_filt_data_neu);
                start_resp = sorted_filt_data_neu(Start_ind,n_row_sorted_filt_data_neu);
                % If the current interval is descending, update the current
                % point to endpoint and update current index to end index
                if curr_resp < start_resp && i_col_sorted_filt_data_neu~=n_col_sorted_filt_data_neu
                    End_point = sorted_filt_data_neu(i_col_sorted_filt_data_neu,1);
                    End_ind = i_col_sorted_filt_data_neu;
                % If the current interval is not descending (or anymore), 
                % compare the current interval with existed maximum
                % interval.
                else
                    % Calculate the current interval
                    curr_interval = End_point-Start_point;
                    % Compare it with the current existed maximum interval
                    % If it's larger, update the existed maximum interval
                    % value.
                    if curr_interval > max_interval || (curr_resp < start_resp && i_col_sorted_filt_data_neu==n_col_sorted_filt_data_neu)
                        max_interval = curr_interval;
                        max_start_point = Start_point;
                        max_end_point = End_point;
                        if i_col_sorted_filt_data_neu==n_col_sorted_filt_data_neu
                            disp('it passed here')
                            guess_x0 = sorted_filt_data_neu(i_col_sorted_filt_data_neu,1);
                        else
                            guess_x0 = (max_start_point+max_end_point)/2;
                        end
                    end
                    Start_point = sorted_filt_data_neu(i_col_sorted_filt_data_neu,1);
                    Start_ind = i_col_sorted_filt_data_neu;
                    End_point = sorted_filt_data_neu(i_col_sorted_filt_data_neu,1);
                    End_ind = i_col_sorted_filt_data_neu;
                end
            end
        end
        lumin_high = max(filt_unique_resp_neu);
        lumin_low = min(filt_unique_resp_neu);
        % Randomly Sample a Point in Data (But Avoid Choosing lumin_high
        % and lumin_low)
        lumin_for_sample = filt_unique_resp_neu(filt_unique_resp_neu~=lumin_high & filt_unique_resp_neu~=lumin_low);
        lumin_for_sample = sort(lumin_for_sample);
        length_lumin_sample = size(lumin_for_sample,1);
        lumin_samp = lumin_for_sample(ceil(length_lumin_sample/2));
        [sam_ind_x,sam_ind_y] = find(filt_unique_resp_neu==lumin_samp);
        conc_sam = filt_unique_pred_neu(sam_ind_x,sam_ind_y);
        if mean_current_NEU_bla>lumin_low
            % Upper Bound Setting
            L = lumin_high-mean_current_NEU_bla; % since the upper bound is supposed to be 100%  
            cal_var1 = log((lumin_high-mean_current_NEU_bla)/(lumin_samp-mean_current_NEU_bla)-1);
        else
            % Upper Bound Setting
            L = lumin_high-lumin_low; % since the upper bound is supposed to be 100%  
            cal_var1 = log((lumin_high-lumin_low)/(lumin_samp-lumin_low)-1);
        end
        cal_var2 = (-1)/(conc_sam-guess_x0);
        k = cal_var1*cal_var2;
        % Determine the minimum (lower) boundary of the function
        if mean_current_NEU_bla>lumin_low
            model3 = @(phi,c)(L./(1+exp(-k.*(c-phi)))+mean_current_NEU_bla);
        else
            model3 = @(phi,c)(L./(1+exp(-k.*(c-phi)))+lumin_low);
        end
        % initial guess
        phi0 = guess_x0; 
        [phi_neu,res_neu] = nlinfit(pred_neu,resp_neu,model3,phi0);
        % Plotting Outcomes
        min_pre = min(pred_neu);
        max_pre = max(pred_neu);
        sim_input = min_pre:0.00001:max_pre;
        % Plotting the model
        figure(pic_counter_NEU + 4*(i_drug-1));
        hold on
        plot(sim_input,model3(phi_neu,sim_input),'k','LineWidth',2)
        out_neu = model3(phi_neu,sim_input);
        % plot display setting
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        hold on
        % Outputing Mean Squared Error
        num0bs = length(pred_neu);
        numParams = 1;
        df = num0bs - numParams;
        mse = (res_neu'*res_neu)/df;
        disp('The mean squared error of the first model on Neuro is:')
        disp(mse)

        % Divided by Area Portion
        % Same idea from HepG2 part: obtained overall neuro data, then
        % divided its luminscence data with its area data
        data_current_NEU_zero = data_current(data_current.conc==0&...
                                             ~contains(data_current.cell,'blank')&...
                                             ~contains(data_current.cell,'-inc')&...
                                              contains(data_current.cell,'neuro'),:);
        data_current_NEU_zero_lumin = data_current_NEU_zero.lumin;
        data_current_NEU_zero_area = data_current_NEU_zero.area;
        data_current_NEU_zero_lumin_divided = data_current_NEU_zero_lumin./data_current_NEU_zero_area;
        % Same idea: calculate the mean
        mean_current_NEU_zero_divided = mean(data_current_NEU_zero_lumin_divided);
        % Same idea: obtains the current celltype area data, then calculate
        % its mean, which will be used as later division of blank
        % luminscence
        data_current_NEU_lumin = data_current.area(contains(data_current.cell,'neuro'),:);
        mean_area_current_NEU = mean(data_current_NEU_lumin);
        data_current_NEU_bla_lumin_divided = data_current_NEU_bla_lumin/mean_area_current_NEU;
        mean_current_NEU_bla_divided = mean(data_current_NEU_bla_lumin_divided);
        mean_current_NEU_bla_divided = mean_current_NEU_bla_divided/mean_current_NEU_zero_divided;
        ci_95_current_NEU_bla_divided = Tbla_multiplier*std(data_current_NEU_bla_lumin_divided)/sqrt(n_bla);
        ci_95_current_NEU_bla_divided = ci_95_current_NEU_bla_divided/mean_current_Hep_zero_divided;
        % Divided by Area Plot
        figure(pic_counter_NEU + 4*(i_drug-1)+1);
        hold on
        upp = (mean_current_NEU_bla_divided  + ci_95_current_NEU_bla_divided);
        yline(upp, '--')
        hold on
        dow = (mean_current_NEU_bla_divided  - ci_95_current_NEU_bla_divided);
        yline(dow, '--')
        hold on
        

        % Model Fitting 
        % Construct Models for Fitting
        % Setting Up Fixed-effect Coefficients
        % Scenario 1: x0 only for random effect

        % Preprocessing the Data
        % Get All Unique Concentration Data
        unique_pred_neu_divided = unique(pred_neu_divided);
        % Create New Storage Variable for Storing Filtered Data
        filt_unique_pred_neu_divided = [];
        filt_unique_resp_neu_divided = [];
        % Looping through Each Unique Concentration
        num_uniq_pred_neu_divided = size(unique_pred_neu_divided,1);
        for i_uniq_pred_neu_divided=1:num_uniq_pred_neu_divided
            curr_unique_pred_neu_divided = unique_pred_neu_divided(i_uniq_pred_neu_divided);
            % Concatenate the current concentration into the new storage
            % variable for predictor
            filt_unique_pred_neu_divided = cat(1,filt_unique_pred_neu_divided,curr_unique_pred_neu_divided);
            % Grab out all indices that are corresponding to this
            % concentration 
            curr_unique_pred_neu_divided_ind = find(pred_neu_divided==curr_unique_pred_neu_divided);
            % Get the frequency of the corresponding concentration
            num_curr_unique_pred_neu_divided = size(curr_unique_pred_neu_divided_ind,1);
            % If there's more than one, average through all the
            % corresponding data and store it into the storage variable
            if num_curr_unique_pred_neu_divided > 1
               avg_curr_unique_neu_resp_divided = mean(resp_neu_divided(curr_unique_pred_neu_divided_ind));
               filt_unique_resp_neu_divided = cat(1,filt_unique_resp_neu_divided, avg_curr_unique_neu_resp_divided);
            % If the corresponding concentration only have one
            % corresponding response, just store it directly.
            else
               filt_unique_resp_neu_divided = cat(1,filt_unique_resp_neu_divided, resp_neu(curr_unique_pred_neu_divided_ind));
            end
        end

        % Looping through each of the elements to check the longest
        % decreasing interval
        % Concatenate Predictor and Response Array Together to Form
        % Corresponding Predictor-Response Pairs Data
        filt_data_neu_divided = cat(2, filt_unique_pred_neu_divided, filt_unique_resp_neu_divided);
        % Sort Out Data Array in Ascending Order of the First Column (Predictor Array)
        % and Grab Out Its Size
        sorted_filt_data_neu_divided = sortrows(filt_data_neu_divided);
        n_col_sorted_filt_data_neu_divided = size(sorted_filt_data_neu_divided,1);
        n_row_sorted_filt_data_neu_divided = size(sorted_filt_data_neu_divided,2);
        % Create Variable for Storing Maximum Interval of Response
        % Descending and Store the Corresponding Start and End predictors 
        max_interval = 0;
        max_start_point = 0;
        max_end_point = 0;
        % Create Variable to Grab Out the Starting and Ending Point of the Current Interval
        Start_point = 0;
        End_point = 0;
        % Create Variable for Storing Starting and Ending Index (for
        % pointing purpose at the data matrix)
        Start_ind = 0;
        End_ind = 0;
        % Loop through the data
        for i_col_sorted_filt_data_neu_divided = 1:n_col_sorted_filt_data_neu_divided
            % Check if the interval is descending
            % If this is the first looping, save the corresponding
            % predictor value and index and skip the investigation
            if i_col_sorted_filt_data_neu_divided == 1
                Start_point = sorted_filt_data_neu_divided(i_col_sorted_filt_data_neu_divided,1);
                Start_ind = i_col_sorted_filt_data_neu_divided;
            % If this is not the first looping, check if the interval is
            % descending
            else
                % Grab out current tested response and the starting point
                curr_resp = sorted_filt_data_neu_divided(i_col_sorted_filt_data_neu_divided,n_row_sorted_filt_data_neu_divided);
                start_resp = sorted_filt_data_neu_divided(Start_ind,n_row_sorted_filt_data_neu_divided);
                % If the current interval is descending, update the current
                % point to endpoint and update current index to end index
                if curr_resp < start_resp && i_col_sorted_filt_data_neu_divided~=n_col_sorted_filt_data_neu_divided
                    End_point = sorted_filt_data_neu_divided(i_col_sorted_filt_data_neu_divided,1);
                    End_ind = i_col_sorted_filt_data_neu_divided;
                % If the current interval is not descending (or anymore), 
                % compare the current interval with existed maximum
                % interval.
                else
                    % Calculate the current interval
                    curr_interval = End_point-Start_point;
                    % Compare it with the current existed maximum interval
                    % If it's larger, update the existed maximum interval
                    % value.
                    if curr_interval > max_interval || (curr_resp < start_resp && i_col_sorted_filt_data_neu_divided==n_col_sorted_filt_data_neu_divided) 
                        max_interval = curr_interval;
                        max_start_point = Start_point;
                        max_end_point = End_point;
                        if i_col_sorted_filt_data_neu_divided==n_col_sorted_filt_data_neu_divided
                            disp('it passed here')
                            guess_x0 = sorted_filt_data_neu_divided(i_col_sorted_filt_data_neu_divided,1);
                        else
                            guess_x0 = (max_start_point+max_end_point)/2;
                        end
                    end
                    Start_point = sorted_filt_data_neu_divided(i_col_sorted_filt_data_neu_divided,1);
                    Start_ind = i_col_sorted_filt_data_neu_divided;
                    End_point = sorted_filt_data_neu_divided(i_col_sorted_filt_data_neu_divided,1);
                    End_ind = i_col_sorted_filt_data_neu_divided;
                end
            end
        end
        lumin_high = max(filt_unique_resp_neu_divided);
        lumin_low = min(filt_unique_resp_neu_divided);
        % Randomly Sample a Point in Data (But Avoid Choosing lumin_high
        % and lumin_low)
        lumin_for_sample = filt_unique_resp_neu_divided(filt_unique_resp_neu_divided~=lumin_high & filt_unique_resp_neu_divided~=lumin_low);
        lumin_for_sample = sort(lumin_for_sample);
        length_lumin_sample = size(lumin_for_sample,1);
        lumin_samp = lumin_for_sample(ceil(length_lumin_sample/2));
        [sam_ind_x,sam_ind_y] = find(filt_unique_resp_neu_divided==lumin_samp);
        conc_sam = filt_unique_pred_neu_divided(sam_ind_x,sam_ind_y);
        if mean_current_NEU_bla_divided>lumin_low
            % Upper Bound Setting
            L = lumin_high-mean_current_NEU_bla_divided; % since the upper bound is supposed to be 100%  
            cal_var1 = log((lumin_high-mean_current_NEU_bla_divided)/(lumin_samp-mean_current_NEU_bla_divided)-1);
        else
            % Upper Bound Setting
            L = 1-lumin_low; % since the upper bound is supposed to be 100%  
            cal_var1 = log((lumin_high-lumin_low)/(lumin_samp-lumin_low)-1);
        end
        cal_var2 = (-1)/(conc_sam-guess_x0);
        k = cal_var1*cal_var2;
        % Determine the minimum (lower) boundary of the function
        if mean_current_NEU_bla_divided>lumin_low
            model4 = @(phi,c)(L./(1+exp(-k.*(c-phi)))+mean_current_NEU_bla_divided);
        else
            model4 = @(phi,c)(L./(1+exp(-k.*(c-phi)))+lumin_low);
        end
        % initial guess
        phi0 = guess_x0; 
        [phi_neu_divided,res_neu_divided] = nlinfit(pred_neu_divided,resp_neu_divided,model4,phi0);
        % Plotting Outcomes
        min_pre = min(pred_neu_divided);
        max_pre = max(pred_neu_divided);
        sim_input = min_pre:0.00001:max_pre;
        % Plotting the model
        figure(pic_counter_NEU + 4*(i_drug-1)+1);
        hold on
        plot(sim_input,model4(phi_neu_divided,sim_input),'k','LineWidth',2)
        out_neu_divided = model4(phi_neu_divided,sim_input);
        % plot display setting
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        hold on
        % Outputing Mean Squared Error
        num0bs = length(pred_neu_divided);
        numParams = 1;
        df = num0bs - numParams;
        mse = (res_neu_divided'*res_neu_divided)/df;
        disp('The mean squared error of the first model on Neuro (Divided) is:')
        disp(mse)
        
        
    end 
         
end