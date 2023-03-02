function up_img=intep_up(img,factor)
    % get the size of the image
    wid_img=size(img,1); % get the width of the image
    len_img=size(img,2); % get the length of the image
    chn_img=size(img,3); % get the channel number of the image
    
    % conducting upsampling
    % conduct upsampling for row direction
    new_img_temp=[]; % create a temporarily array for storage
    for i_wid=1:wid_img
        % condition which channel is not single
        % extract each row of the image
        if chn_img ~= 1
            % create a temporarily array for storing row data
            new_img_temp_row=[];
            % looping through each channel
            for i_chn=1:chn_img
                % grab out data for each row in each channel
                curr_img=img(i_wid,:,i_chn);
                % conducting upsampling
                new_curr_img=interp(double(curr_img),factor);
                % concatenate the data along channel direction
                new_img_temp_row=cat(3,new_img_temp_row,new_curr_img);
            end
            % concatenate the row data for all channel to temporarily
            % storage for the entire image
            new_img_temp=[new_img_temp; new_img_temp_row];
        
        % single channel image case
        else
            % grab out row data for the image
            curr_img=img(i_wid,:);
            % conducting upsampling
            new_curr_img=interp(double(curr_img),facotr);
            new_img_temp=[new_img_temp; new_curr_img];
        end
    end

    % conduct upsampling for column direction
    new_img=[];
    n_col=size(new_img_temp,2);
    % looping through each column
    for i_col=1:n_col
        % condution which channel is not single
        % extract each column of the image
        if chn_img ~= 1
            % create a temporarily array for storing column data
            new_img_temp_col=[];
            % looping through each channel
            for i_chn=1:chn_img
                % grab out data for each column in each channel
                curr_img=new_img_temp(:,i_col,i_chn);
                % conducting upsampling 
                new_curr_img=interp(double(curr_img),factor);
                % concatenate the data along channel direction
                new_img_temp_col=cat(3,new_img_temp_col,new_curr_img);
            end
            % catenate the column data for all channel to temporarily
            % storage for the entire image
            new_img=[new_img new_img_temp_col];

        % single channel image case
        else
            % grab out column data for the image
            curr_img=new_img_temp(:,i_col);
            % conduct upsampling
            new_curr_img=interp(double(curr_img),factor);
            new_img=[new_img new_curr_img];
        end
    end
    up_img=uint8(new_img);
end