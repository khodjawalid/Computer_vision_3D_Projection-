function frames = extractFrames(videoFile)
%This function allows us to extract frames from a video stream

% Load the video
vid=VideoReader(videoFile);

% read the total number of frames
nbFrames = vid.NumberOfFrames;

% fps = get(videoReader, 'FrameRate');
% disp(fps); % the fps is correct: it's the same declared in the video file properties

frames = cell(vid.NumFrames,1);
i = 1;
while hasFrame(vid)
    frames{i} = readFrame(vid);
    i = i + 1;
end 


% saving the extracted frames in a folder 
Folder = 'myVideoFrames2/';
for i = 1:nbFrames
    iFrame = read(vid, i);
    imwrite(iFrame, fullfile(Folder, sprintf('%06d.jpg', i)));
end 

end


