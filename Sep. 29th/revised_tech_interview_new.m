function revised_tech_interview_new(input_filename, input_celltypes)
    
    path(path,'../../matlab-lib');
    
    %% Check Out All Input 
    if isempty(input_filename)
        error("No input filename found. Please enter the intended analyzed file")
    end
    
    
    
    %% Input Data
    % Import Filename 
    filename = input_filename;
    % Find Out How Many Spreadsheets
    [~, sheets] = xlsfinfo(filename);
        
    % Parameters Settings For Further Purposes
    % For Looping Purpose
    n_drugs = length(sheets);
    n_cells = numel(input_celltypes);
    
    % Determine how many analysis you would like to perform
    analy_type = ["Original" "Divided by Area" "Divided by Volume"];
    num_analy = numel(analy_type);
    
%      ## this type of hard coding is not good

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
        data_current.uid = string(data_current.uid);
        data_current.did = string(data_current.did);
        data_current.drug = string(data_current.drug);
        data_current.cell = string(data_current.cell);
        
        % Check If All the Input Celltypes Are Inside the Dataset
        % Since the celltype is an input variable, there might be a typo. 
        % This section can check if the cell type indicated in the input is
        % also within the data
        for i_inp_cell=1:length(input_celltypes)
           if ~ismember(string(input_celltypes(i_inp_cell)),data_current.cell)
               message = "Uknown cellType found : %s. Existed cellType: %s";
               error(message,string(input_celltypes(i_inp_cell)),data_current.cell);
           end
        end
        
%         ## it's a good practice to keep all string information of a table in string type.
        % Create Variable for Storing Predictors and Responses
%         ## the below table data_in_ids has both HepG2 and neuro, so the line color is alternatively used.  How to fix this?
        for i_cell=1:n_cells
            % Grab out current cell type
            celltype = input_celltypes(i_cell);
            % Grab out only the data for current cell type
            data_current_cell = data_current(contains(data_current.cell, celltype), :);
%             ## filter the table so that it has only 1st cell type, what ever it is
%             
%             ## we reset these variable when we start a new celltype
            
            % Setup Portion
            % Storage Setup
            % Create Corresponding Varibles for Corresponding Analysis
            % Method
            % For current analysis, it contains for 
            % 1. Original Data
            % 2. Data divided Area
            % 3. Data divided Volume
            for i_method=1:num_analy
                research_method = analy_type(i_method);
                % do nothing. use the original value
                if strcmp(research_method,"Original")
                    pred = [];
                    resp = [];
                % normalize by the area
                elseif strcmp(research_method,"Divided by Area")
                    pred_divided = [];
                    resp_divided = [];
                % normalize by the volume
                elseif strcmp(research_method,"Divided by Volume")
                    pred_divided_vol = [];
                    resp_divided_vol = [];
                % print out an error message is the input analysis method
                % is unknown
                else
                    message = "New Research Method %s Has Been Added.";
                    method = research_method;
                    error(message,method); 
                end
                
                % Parameters Setup
                % Calculate the Number of Data IDs in This SpreadSheet
                % Number of Data IDs Calculation
                % Retrieve the data id for each sheet and calculate using numel
                % function
                data_current_cell_ids = unique(data_current_cell.did);
                data_current_cell_num_dids = numel(data_current_cell_ids); 
                
                % Plot Setup
                pic_counter = (i_cell-1)*num_analy+i_method;
                figure(pic_counter);
                message = "Dose vs Concentration Plot for %s (%s)";
                title(sprintf(message, celltype, research_method));
                clr = lines(data_current_cell_num_dids);
                % Plot display setting
                set(gca, "XScale", "log")
                set(gca, "YScale", "log")
                % Label the plot
                xlabel("Concentration")
                ylabel("Luminance (%)")
                grid on
                hold on
                % Setup Fake X Point for Zero Concentration
                % Set up a fake x axis for zero concentration by finding the
                % minimum value of all nonzero concentration in the spread sheet and divided
                % it by 10^6
                conc_min = min(data_current_cell.conc(data_current_cell.conc>0,:));
                factor = 1e6;
                fake_x_zero = conc_min/factor;
                % Setup Fake X Point for Incubator
                % Set up a fake x axis for zero concentration by finding the
                % minimum value of all nonzero concentration in the spread sheet and divided
                % it by 10^7            
                factor = 1e7;
                fake_x_inc = conc_min/factor;
                
                % Data Filtering Portion
                % Looping through Each Batch 
                for i_did=1:data_current_cell_num_dids
                   % Filter Out Unnecessary Data 
                   % (preserve only the data in current batch)
                   did_current_cell_did = data_current_cell_ids(i_did);
                   data_current_cell_did = data_current_cell(data_current_cell.did == did_current_cell_did,:);
                   % Data Retreiving Portion 
                   % Zero Concentration Data
                   % (For filtering out data that are blank or incubator cases)
                   data_current_cell_did_zero = data_current_cell_did(data_current_cell_did.conc == 0 & ...
                                                                      ~contains(data_current_cell_did.cell,"-inc"),:);
                   data_current_cell_did_zero_lumin = data_current_cell_did_zero.lumin;
                   % Incubator Data
                   % (For filtering out data that are zero concentration or blank cases)
                   data_current_cell_did_inc = data_current_cell_did(contains(data_current_cell_did.cell, "-inc"),:);
                   data_current_cell_did_inc_lumin = data_current_cell_did_inc.lumin;                 

                   % Calculation Portion
                   % Original Portion
                   if strcmp(research_method,"Original")
                        % (Make sure it's divided by itself to make it percentile)
                        % Calculate the Mean and CI_95 for Zero Concentration Data
                        mean_current_cell_did_zero_temp = mean(data_current_cell_did_zero_lumin);
                        ci_95_current_cell_did_zero = SEtoCI(std(data_current_cell_did_zero_lumin), numel(data_current_cell_did_zero_lumin), true);
                        % normalize mean and ci 
                        mean_current_cell_did_zero = mean_current_cell_did_zero_temp/mean_current_cell_did_zero_temp;
                        ci_95_current_cell_did_zero = ci_95_current_cell_did_zero/mean_current_cell_did_zero_temp;
                        % Calculate the Mean and CI_95 for Incubator Data
                        mean_current_cell_did_inc = mean(data_current_cell_did_inc_lumin);
                        ci_95_current_cell_did_inc = SEtoCI(std(data_current_cell_did_inc_lumin), numel(data_current_cell_did_inc_lumin), true);
                        % normalize mean and ci 
                        mean_current_cell_did_inc = mean_current_cell_did_inc/mean_current_cell_did_zero_temp;
                        ci_95_current_cell_did_inc = ci_95_current_cell_did_inc/mean_current_cell_did_zero_temp;

                        % Plotting Portion
                        % Plotting Both Zero and Incubator in the Current Plot
                        % Plotting Zero Concentration
                        e = errorbar(fake_x_zero,mean_current_cell_did_zero,ci_95_current_cell_did_zero,"color",clr(i_did,:), "linestyle","none");
                        e.LineStyle = "none";
                        hold on
                        % Plotting Incubaotr Concentration
                        e = errorbar(fake_x_inc,mean_current_cell_did_inc,ci_95_current_cell_did_inc,"color",clr(i_did,:), "linestyle","none");
                        e.LineStyle = "none";
                        hold on

                        % Storing Portion
                        % Store Data for Later Model Fitting Purpose
                        if isempty(pred) && isempty(resp)
                          pred = cat(1,fake_x_zero,fake_x_inc);
                          resp = cat(1,mean_current_cell_did_zero,mean_current_cell_did_inc);
                        elseif ~isempty(pred) && ~isempty(resp)
                          pred = cat(1,pred,fake_x_zero,fake_x_inc);
                          resp = cat(1,resp,mean_current_cell_did_zero,mean_current_cell_did_inc);
                        else
                           % ## HW: same here, use sprintf()
                           message = "%s (%s) Goes Wrong. Unmatched Data Storage";
                           error(message, celltype, research_method);
                        end

                    % Divided by Area Portion
                    elseif strcmp(research_method,"Divided by Area")
                        % (Make sure it's divided by itself to make it percentile)
                        % Retrieve Area Data
                        data_current_cell_did_zero_area = data_current_cell_did_zero.area;
                        % Divide the Data with Area
                        data_current_cell_did_zero_lumin_divided = data_current_cell_did_zero_lumin./data_current_cell_did_zero_area;
                        % Calculate the Mean and CI_95 for Zero Concentration Data
                        mean_current_cell_did_zero_divided_temp = mean(data_current_cell_did_zero_lumin_divided);
                        ci_95_current_cell_did_zero_divided = SEtoCI(std(data_current_cell_did_zero_lumin_divided), numel(data_current_cell_did_zero_lumin_divided), true);
                        % normalize mean and ci 
                        mean_current_cell_did_zero_divided = mean_current_cell_did_zero_divided_temp/mean_current_cell_did_zero_divided_temp;
                        ci_95_current_cell_did_zero_divided = ci_95_current_cell_did_zero_divided/mean_current_cell_did_zero_divided_temp;
                        % Retrieve Area Data
                        data_current_cell_did_inc_area = data_current_cell_did_inc.area;
                        % Divide the Data with Area
                        data_current_cell_did_inc_lumin_divided = data_current_cell_did_inc_lumin./data_current_cell_did_inc_area;
                        % Calculate the Mean and CI_95 for Incubator Data
                        mean_current_cell_did_inc_divided = mean(data_current_cell_did_inc_lumin_divided);
                        ci_95_current_cell_did_inc_divided = SEtoCI(std(data_current_cell_did_inc_lumin_divided), numel(data_current_cell_did_inc_lumin_divided), true);
                        % normalize mean and ci
                        mean_current_cell_did_inc_divided = mean_current_cell_did_inc_divided/mean_current_cell_did_zero_divided_temp;
                        ci_95_current_cell_did_inc_divided = ci_95_current_cell_did_inc_divided/mean_current_cell_did_zero_divided_temp;

                        % Plotting Portion
                        % Plotting Both Zero and Incubator in the Current Plot
                        % Plotting Zero Concentration
                        e = errorbar(fake_x_zero,mean_current_cell_did_zero_divided,ci_95_current_cell_did_zero_divided,"color",clr(i_did,:), "linestyle","none");
                        e.LineStyle = "none";
                        hold on 
                        % Plotting Incubaotr Concentration
                        e = errorbar(fake_x_inc,mean_current_cell_did_inc_divided,ci_95_current_cell_did_inc_divided,"color",clr(i_did,:), "linestyle","none");
                        e.LineStyle = "none";
                        hold on

                        % Storing Portion
                        if isempty(pred_divided) && isempty(resp_divided)
                           pred_divided = cat(1,fake_x_zero,fake_x_inc);
                           resp_divided = cat(1,mean_current_cell_did_zero_divided,mean_current_cell_did_inc_divided);
                        elseif ~isempty(pred_divided) && ~isempty(resp_divided)
                           pred_divided = cat(1,pred_divided,fake_x_zero,fake_x_inc);
                           resp_divided = cat(1,resp_divided,mean_current_cell_did_zero_divided,mean_current_cell_did_inc_divided);
                        else
                           % ## HW: same here, use sprintf()
                           message = "%s (%s) Goes Wrong. Unmatched Data Storage";
                           error(message, celltype, research_method);
                        end

                    % Divided by Volume Portion
                    elseif strcmp(research_method,"Divided by Volume") 
                         % (Make sure it's divided by itself to make it percentile)
                         % Calculate the Volume Data
                         % Retrieve Area Data
                         data_current_cell_did_zero_area = data_current_cell_did_zero.area;
                         data_current_cell_did_zero_vol = data_current_cell_did_zero_area./pi;
                         data_current_cell_did_zero_vol = sqrt(data_current_cell_did_zero_vol);
                         data_current_cell_did_zero_vol = power(data_current_cell_did_zero_vol,3);
                         data_current_cell_did_zero_vol = (4/3).*data_current_cell_did_zero_vol;
                         % Divide the Data with Volume
                         data_current_cell_did_zero_lumin_vol = data_current_cell_did_zero_lumin./data_current_cell_did_zero_vol;
                         % Calculate the Mean and CI_95 for Zero Concentration Data
                         mean_current_cell_did_zero_vol_temp = mean(data_current_cell_did_zero_lumin_vol);
                         ci_95_current_cell_did_zero_vol = SEtoCI(std(data_current_cell_did_zero_lumin_vol), numel(data_current_cell_did_zero_lumin_vol), true);
                         % normalize mean and ci 
                         mean_current_cell_did_zero_vol = mean_current_cell_did_zero_vol_temp/mean_current_cell_did_zero_vol_temp;
                         ci_95_current_cell_did_zero_vol = ci_95_current_cell_did_zero_vol/mean_current_cell_did_zero_vol_temp;
                         % Calculate the Volume Data
                         % Retrieve Area Data
                         data_current_cell_did_inc_area = data_current_cell_did_inc.area;
                         data_current_cell_did_inc_vol = data_current_cell_did_inc_area./pi;
                         data_current_cell_did_inc_vol = sqrt(data_current_cell_did_inc_vol);
                         data_current_cell_did_inc_vol = power(data_current_cell_did_inc_vol,3);
                         data_current_cell_did_inc_vol = (4/3).*data_current_cell_did_inc_vol;
                         % Divide the Data with Area
                         data_current_cell_did_inc_lumin_vol = data_current_cell_did_inc_lumin./data_current_cell_did_inc_vol;
                         % Calculate the Mean and CI_95 for Incubator Data
                         mean_current_cell_did_inc_vol = mean(data_current_cell_did_inc_lumin_vol);
                         ci_95_current_cell_did_inc_vol = SEtoCI(std(data_current_cell_did_inc_lumin_vol), numel(data_current_cell_did_inc_lumin_vol), true);
                         % normalize mean and ci
                         mean_current_cell_did_inc_vol = mean_current_cell_did_inc_vol/mean_current_cell_did_zero_vol_temp;
                         ci_95_current_cell_did_inc_vol = ci_95_current_cell_did_inc_vol/mean_current_cell_did_zero_vol_temp;

                         % Plotting Portion
                         % Plotting Both Zero and Incubator in the Current Plot
                         % Plotting Zero Concentration
                         e = errorbar(fake_x_zero,mean_current_cell_did_zero_vol,ci_95_current_cell_did_zero_vol,"color",clr(i_did,:), "linestyle","none");
                         e.LineStyle = "none";
                         hold on 
                         % Plotting Incubaotr Concentration
                         e = errorbar(fake_x_inc,mean_current_cell_did_inc_vol,ci_95_current_cell_did_inc_vol,"color",clr(i_did,:), "linestyle","none");
                         e.LineStyle = "none";
                         hold on

                         % Storing Portion
                         if isempty(pred_divided_vol) && isempty(resp_divided_vol)
                            pred_divided_vol = cat(1,fake_x_zero,fake_x_inc);
                            resp_divided_vol = cat(1,mean_current_cell_did_zero_vol,mean_current_cell_did_inc_vol);
                         elseif ~isempty(pred_divided_vol) && ~isempty(resp_divided_vol)
                            pred_divided = cat(1,pred_divided_vol,fake_x_zero,fake_x_inc);
                            resp_divided = cat(1,resp_divided_vol,mean_current_cell_did_zero_vol,mean_current_cell_did_inc_vol);
                         else
                            % ## HW: same here, use sprintf()
                            message = "%s (%s) Goes Wrong. Unmatched Data Storage";
                            error(message, celltype, research_method);
                         end

                   % If the method indicated is not in the current function,
                   % print out an error message
                   else
                        message = "Research method %s not found. Function needs modification.";
                        error(message,research_method);
                   end


                   % Data Filtering
                   % Grab Out Only Nonzero Concentration Data 
                   data_current_cell_did_concs = data_current_cell_did(data_current_cell_did.conc>0,:);
                   % Detect How Many Concentrations Were Tested in This Batch
                   data_current_cell_did_concs_c = unique(data_current_cell_did_concs.conc); 
                   data_current_cell_did_concs_num = numel(data_current_cell_did_concs_c);
                   % Looping through each concentration
                   for i_conc=1:data_current_cell_did_concs_num
                      % Plot Setup
                      % Setup X point for current concentration
                      current_conc = data_current_cell_did_concs_c(i_conc);
                      % Data Retrieving Portion
                      % Grab Out the Current Concentration Data
                      data_current_cell_did_currconc = data_current_cell_did_concs(data_current_cell_did_concs.conc...
                                                                              ==data_current_cell_did_concs_c(i_conc),:);
                      data_current_cell_did_currconc_lumin = data_current_cell_did_currconc.lumin;

                      % Calculation Portion
                      % Original Data Portion
                      if strcmp(research_method,"Original")
                          % (Make sure it's divided by zero concentration data to make it percentile)
                          % Calculate the Mean and CI_95 for Current Concentration Data
                          mean_current_cell_did_currconc = mean(data_current_cell_did_currconc_lumin);
                          ci_95_current_cell_did_currconc = SEtoCI(std(data_current_cell_did_currconc_lumin), numel(data_current_cell_did_currconc_lumin), true);
                          % normalize mean and ci
                          mean_current_cell_did_currconc = mean_current_cell_did_currconc/mean_current_cell_did_zero_temp;
                          ci_95_current_cell_did_currconc = ci_95_current_cell_did_currconc/mean_current_cell_did_zero_temp;

                          % Plotting Portion
                          % Plot Current Concentration Data
                          e = errorbar(current_conc,mean_current_cell_did_currconc,ci_95_current_cell_did_currconc,"color",clr(i_did,:), "linestyle","none");
                          e.LineStyle = "none";
                          hold on

                          % Storing Portion
                          % Store Data for Later Model Fitting Purpose
                          if isempty(pred) && isempty(resp)
                            pred = current_conc;
                            resp = mean_current_cell_did_currconc;
                          elseif ~isempty(pred) && ~isempty(resp)
                            pred = cat(1,pred,current_conc);
                            resp = cat(1,resp,mean_current_cell_did_currconc);
                          else
                            % ## HW: same here, use sprintf()
                            message = "%s (%s) Goes Wrong. Unmatched Data Storage";
                            error(message, celltype, research_method);
                          end

                      % Divided by Area Portion
                      elseif strcmp(research_method,"Divided by Area")
                         % (Make sure it's divided by itself to make it percentile)
                         % Retrieve Area Data
                         data_current_cell_did_currconc_area = data_current_cell_did_currconc.area;
                         % Divide the Data with Area
                         data_current_cell_did_currconc_lumin_divided = data_current_cell_did_currconc_lumin./data_current_cell_did_currconc_area;
                         % Calculate the Mean and CI_95 for Current
                         % Concentration Data
                         mean_current_cell_did_currconc_divided = mean(data_current_cell_did_currconc_lumin_divided);
                         ci_95_current_cell_did_currconc_divided = SEtoCI(std(data_current_cell_did_currconc_lumin_divided), numel(data_current_cell_did_currconc_lumin_divided), true);
                         % normalize mean and ci
                         mean_current_cell_did_currconc_divided = mean_current_cell_did_currconc_divided/mean_current_cell_did_zero_divided_temp;
                         ci_95_current_cell_did_currconc_divided = ci_95_current_cell_did_currconc_divided/mean_current_cell_did_zero_divided_temp;

                         % Plotting Portion
                         e = errorbar(current_conc,mean_current_cell_did_currconc_divided,ci_95_current_cell_did_currconc_divided,"color",clr(i_did,:), "linestyle","none");
                         e.LineStyle = "none";
                         hold on

                        % Storing Portion
                        if isempty(pred_divided) && isempty(resp_divided)
                           pred_divided = current_conc;
                           resp_divided = mean_current_cell_did_currconc_divided;
                        elseif ~isempty(pred_divided) && ~isempty(resp_divided)
                           pred_divided = cat(1,pred_divided,current_conc);
                           resp_divided = cat(1,resp_divided,mean_current_cell_did_currconc_divided);
                        else
                           % ## HW: same here, use sprintf()
                           message = "%s (%s) Goes Wrong. Unmatched Data Storage";
                           error(message, celltype, research_method);
                        end

                     % Divided by Volume Portion
                     elseif strcmp(research_method,"Divided by Volume") 
                         % (Make sure it's divided by itself to make it percentile)
                         % Calculate the Volume Data
                         % Retrieve Area Data
                         data_current_cell_did_currconc_area = data_current_cell_did_currconc.area;
                         data_current_cell_did_currconc_vol = data_current_cell_did_currconc_area./pi;
                         data_current_cell_did_currconc_vol = sqrt(data_current_cell_did_currconc_vol);
                         data_current_cell_did_currconc_vol = power(data_current_cell_did_currconc_vol,3);
                         data_current_cell_did_currconc_vol = (4/3).*data_current_cell_did_currconc_vol;
                         % Divide the Data with Area
                         data_current_cell_did_currconc_lumin_vol = data_current_cell_did_currconc_lumin./data_current_cell_did_currconc_vol;
                         % Calculate the Mean and CI_95 for Incubator Data
                         mean_current_cell_did_currconc_vol = mean(data_current_cell_did_currconc_lumin_vol);
                         ci_95_current_cell_did_currconc_vol = SEtoCI(std(data_current_cell_did_currconc_lumin_vol), numel(data_current_cell_did_currconc_lumin_vol), true);
                         % normalize mean and ci
                         mean_current_cell_did_currconc_vol = mean_current_cell_did_currconc_vol/mean_current_cell_did_zero_vol_temp;
                         ci_95_current_cell_did_currconc_vol = ci_95_current_cell_did_currconc_vol/mean_current_cell_did_zero_vol_temp;

                         % Plotting Portion
                         e = errorbar(current_conc,mean_current_cell_did_currconc_vol,ci_95_current_cell_did_currconc_vol,"color",clr(i_did,:), "linestyle","none");
                         e.LineStyle = "none";
                         hold on

                         % Storing Portion
                         if isempty(pred_divided_vol) && isempty(resp_divided_vol)
                            pred_divided_vol = current_conc;
                            resp_divided_vol = mean_current_cell_did_currconc_vol;
                         elseif ~isempty(pred_divided_vol) && ~isempty(resp_divided_vol)
                            pred_divided_vol = cat(1,pred_divided_vol,current_conc);
                            resp_divided_vol = cat(1,resp_divided_vol,mean_current_cell_did_currconc_vol);
                         else
                            % ## HW: same here, use sprintf()
                            message = "%s (%s) Goes Wrong. Unmatched Data Storage";
                            error(message, celltype, research_method);
                         end
                      % If the method indicated is not in the current function,
                      % print out an error message
                      else
                        message = "Research method %s not found. Function needs modification.";
                        error(message,research_method);
                      end
                   end      
                end

                % Data Filtering
                % (Make sure it's divided by the overall zero concentration mean)
                data_current_cell_zero = data_current_cell(data_current_cell.cell==string(celltype) &...
                                                           data_current_cell.conc==0,:);
                data_current_cell_zero_lumin = data_current_cell_zero.lumin;
                mean_current_cell_zero = mean(data_current_cell_zero_lumin);
                % Grab Out All Blank Data
                % Create an Indicator for Filrering Out Blank Data
                if strcmp(celltype,"HepG2")
                    indicator = "HEP";
                elseif strcmp(celltype, "neuro")
                    indicator = "NEU";
                else
                    message = "%s is not found";
                    error(message,celltype);
                end
                data_current_cell_bla = data_current(strcmp(data_current.cell,"blank") &...
                                                     contains(data_current.did,indicator),:);
                data_current_cell_bla_lumin = data_current_cell_bla.lumin;

                % Calculation Portion
                % Original Data Portion
                if strcmp(research_method,"Original")
                   % (Make sure it's divided by zero concentration data to make it percentile)
                   % Calculate the Mean and CI_95 for Current Concentration Data
                   mean_current_cell_bla = mean(data_current_cell_bla_lumin);
                   ci_95_current_cell_bla = SEtoCI(std(data_current_cell_bla_lumin), numel(data_current_cell_bla_lumin), true);
                   % normalize mean and ci
                   mean_current_cell_bla = mean_current_cell_bla/mean_current_cell_zero;
                   ci_95_current_cell_bla = ci_95_current_cell_bla/mean_current_cell_zero;

                   % Plotting Portion
                   upp = (mean_current_cell_bla + ci_95_current_cell_bla);
                   yline(upp, '--')
                   hold on
                   dow = (mean_current_cell_bla - ci_95_current_cell_bla);
                   yline(dow, '--')
                   hold on

                   % Model Fitting 
                   % Scenario 1: x0 only for random effect
                   [phi_cell mse model] = modelfitting1(pred,resp,mean_current_cell_bla,celltype);
                   % Plotting Outcomes
                   min_pre = min(pred);
                   max_pre = max(pred);
                   sim_input = min_pre:0.00001:max_pre;
                   % Plotting the model
                   plot(sim_input,model(phi_cell,sim_input),'k','LineWidth',2)
                   hold on

                % Divided by Area Portion
                elseif strcmp(research_method,"Divided by Area")
                    % Calculation Portion
                    % In addition to the area of current celltype, the mean of divided zero
                    % concentration data is also needed
                    % To get that, we first need to grab out all zero concentration in
                    % this celltype and then divided its luminscences by its areas 
                    data_current_cell_zero_area = data_current_cell_zero.area;
                    data_current_cell_zero_divided = data_current_cell_zero_lumin./data_current_cell_zero_area;
                    % After that, we calculate its mean
                    mean_current_cell_zero_divided = mean(data_current_cell_zero_divided);
                    % For blank cases, it needs to be divided by the mean area value of
                    % the current celltype.
                    % We then first grab out the area data and calculate its mean
                    data_current_cell_area = data_current.area(contains(data_current.cell,string(celltype)),:);
                    mean_area_current_cell = mean(nonzeros(data_current_cell_area));
                    % After doing so, we then divided the blank data with the mean area
                    % of this cellType
                    data_current_cell_bla_lumin_divided = data_current_cell_bla_lumin/mean_area_current_cell;
                    mean_current_cell_bla_divided = mean(data_current_cell_bla_lumin_divided);
                    ci_95_current_cell_bla_divided = SEtoCI(std(data_current_cell_bla_lumin_divided), numel(data_current_cell_bla_lumin_divided), true);
                    % To achieve the final percentile data, we then utilize the mean of
                    % overall divided zero concentration and divided our obtained mean blank value by it 
                    mean_current_cell_bla_divided = mean_current_cell_bla_divided/mean_current_cell_zero_divided;
                    ci_95_current_cell_bla_divided = ci_95_current_cell_bla_divided/mean_current_cell_zero_divided;

                    % Plotting Portion
                    upp = (mean_current_cell_bla_divided + ci_95_current_cell_bla_divided);
                    yline(upp, '--')
                    hold on
                    dow = (mean_current_cell_bla_divided - ci_95_current_cell_bla_divided);
                    yline(dow, '--')
                    hold on

                    % Model Fitting 
                    [phi_cell_divided mse_divided model_divided] = modelfitting1(pred_divided, resp_divided, mean_current_cell_bla_divided, celltype);
                    % Plotting Outcomes
                    min_pre = min(pred_divided);
                    max_pre = max(pred_divided);
                    sim_input = min_pre:0.00001:max_pre;
                    % Plotting the model
                    plot(sim_input,model_divided(phi_cell_divided,sim_input),'k','LineWidth',2)
                    hold on

                 % Divided by Volume Portion
                 elseif strcmp(research_method,"Divided by Volume") 
                    % Calculation Portion
                    % In addition to the area of current celltype, the mean of divided volume zero
                    % concentration data is also needed
                    % Calculate the Volume of the Zero Concentration Data
                    data_current_cell_zero_area = data_current_cell_zero.area;
                    data_current_cell_zero_vol = data_current_cell_zero_area./pi;
                    data_current_cell_zero_vol = sqrt(data_current_cell_zero_vol);
                    data_current_cell_zero_vol = power(data_current_cell_zero_vol,3);
                    data_current_cell_zero_vol = (4/3*pi).*data_current_cell_zero_vol;
                    data_current_cell_zero_divided_vol = data_current_cell_zero_lumin./data_current_cell_zero_vol;
                    % After that, we calculate its mean
                    mean_current_cell_zero_divided_vol = mean(data_current_cell_zero_divided_vol);
                    % For blank cases, it needs to be divided by the mean area value of
                    % the current celltype.
                    % We then first grab out the area data 
                    data_current_cell_area = data_current.area(contains(data_current.cell,string(celltype)),:);
                    % Calculate the Corresponding Volume
                    data_current_cell_vol = data_current_cell_area./pi;
                    data_current_cell_vol = sqrt(data_current_cell_vol);
                    data_current_cell_vol = power(data_current_cell_vol,3);
                    data_current_cell_vol = (4/3*pi).* data_current_cell_vol;
                    mean_vol_current_cell = mean(nonzeros(data_current_cell_vol));
                    % After doing so, we then divided the blank data with the
                    % mean volume of this cellType
                    data_current_cell_bla_lumin_divided_vol = data_current_cell_bla_lumin/mean_vol_current_cell;
                    mean_current_cell_bla_divided_vol = mean(data_current_cell_bla_lumin_divided_vol);
                    ci_95_current_cell_bla_divided_vol = SEtoCI(std(data_current_cell_bla_lumin_divided_vol), numel(data_current_cell_bla_lumin_divided_vol), true);
                    % To achieve the final percentile data, we then utilize the mean of
                    % overall divided zero concentration and divided our obtained mean blank value by it 
                    mean_current_cell_bla_divided_vol = mean_current_cell_bla_divided_vol/mean_current_cell_zero_divided_vol;
                    ci_95_current_cell_bla_divided_vol = ci_95_current_cell_bla_divided_vol/mean_current_cell_zero_divided_vol;

                    % Plotting Portion
                    upp = (mean_current_cell_bla_divided_vol + ci_95_current_cell_bla_divided_vol);
                    yline(upp, '--')
                    hold on
                    dow = (mean_current_cell_bla_divided_vol - ci_95_current_cell_bla_divided_vol);
                    yline(dow, '--')
                    hold on

                    % Model Fitting 
                    [phi_cell_divided_vol mse_divided_vol model_divided_vol] = modelfitting1(pred_divided_vol, resp_divided_vol, mean_current_cell_bla_divided_vol, celltype);
                    % Plotting Outcomes
                    min_pre = min(pred_divided);
                    max_pre = max(pred_divided);
                    sim_input = min_pre:0.00001:max_pre;
                    % Plotting the model
                    plot(sim_input,model_divided_vol(phi_cell_divided_vol,sim_input),'k','LineWidth',2)
                    hold on 
                    
                 else
                    message = "Research method %s not found. Function needs modification.";
                    error(message,research_method);
                 end
            % The end for Loop: Analysis Method
            end
        % The end for Loop: Celltype
        end 
    % The end for Loop: Drug Type
    end
    
% The end for the function    
end