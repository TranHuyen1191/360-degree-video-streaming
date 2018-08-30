function [Y, U, V, read_fr_id, lost_fr_count] = F_loadFileYuv(fileName, width, height,total_fr_num, start_fr_id, num_fr_read, lost_fr_list)
fileId = fopen(fileName, 'r');
% disp(total_fr_num);
Y=zeros(height, width,num_fr_read);
U=zeros(height/2, width/2,num_fr_read);
V=zeros(height/2, width/2,num_fr_read);

%
read_fr_id = 1;
fr_id = start_fr_id;
lost_fr_count = 0; 

% seek to the start_fr_id position
frame_size = 1.5 * width * height;
fseek(fileId, (start_fr_id -1) * frame_size, 'bof');

% Run on original frame order
while fr_id <= total_fr_num && read_fr_id <= num_fr_read
   
    
    buf = fread(fileId, width * height, 'uchar');
    
    if ismember(fr_id, lost_fr_list) == 0
        % read Y component
        Y(:, :, read_fr_id) = reshape(buf, width, height).'; % reshape
         % read U component
        buf = fread(fileId, width / 2 * height / 2, 'uchar');
        U(:, :, read_fr_id) = reshape(buf, width/2, height/2).'; % reshape
        % read V component
        buf = fread(fileId, width / 2 * height / 2, 'uchar');
        V(:, :, read_fr_id) = reshape(buf, width/2, height/2).'; % reshape
        
        read_fr_id = read_fr_id + 1;
    else
        lost_fr_count = lost_fr_count + 1;
        fprintf('lost_fr %d: %d\n',lost_fr_count, fr_id);
    end
    
   
    % 
    fr_id = fr_id + 1;

end
read_fr_id = read_fr_id - 1;

fclose(fileId);