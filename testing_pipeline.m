%% Settings

scale=0.8720;
sensitivity_bin=0.51;

%% Choose image from directory

myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.tif')); %gets all wav files in struct

%Read files in and create 3D image stack
for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    fullFileName = fullfile(myDir, baseFileName);
    x = imread(fullFileName, 'tif');
    if k==1
       dimensions=size(x,[1 2]);
       array=zeros([dimensions,size(myFiles)],"uint8");
    end 
    gray= rgb2gray(x);
    array(:,:,k)=gray;
end

% Adapted from:
% https://uk.mathworks.com/matlabcentral/answers/29837-processing-files-using-a-for-loop

%% Select slice and select bead to analyse

% Crop image
slice_selected=14;
imshow(array(:,:,slice_selected))
str={'Slice',slice_selected-1};
title(str)

roi = drawrectangle("FixedAspectRatio",true);
loc= roi.Position;
loc=floor(loc);

imshow(array(:,:,slice_selected))
rectangle('Position',loc)

cropped_array=array(loc(2):(loc(2)+loc(4)),loc(1):(loc(1)+loc(3)),slice_selected);
imshow(cropped_array)

%% Detect the bead and measure its diameter

[centers,radii]=detect_beads(cropped_array); % HARDCODED + radius inside function
figure;
imshow(cropped_array);
h = viscircles(centers,radii);
%diameter=2*(radii/0.88);
if size(centers,1)>0
    diameter=2*(radii/scale); % HARDCODED
    writematrix(diameter,'diameters.csv','WriteMode','append')
end
%clear("centers","radii");
pause(0.05);

if isempty(centers) 
    roi_bead = drawcircle;
    centers=roi_bead.Center;
    radii=roi_bead.Radius;
    diameter=2*(radii/scale);
end

%% Transform Image to Fourier Space

input=double(cropped_array);
j = imflatfield(input,60);
J=fftshift(fft2(j));

% Apply Low Pass Filter
low_pass=zeros(size(cropped_array,1));
low_pass_size=floor(size(cropped_array,1)/4);
low_pass((size(cropped_array,1)/2-low_pass_size):(size(cropped_array,1)/2+low_pass_size),(size(cropped_array,1)/2-low_pass_size):(size(cropped_array,1)/2+low_pass_size)) = 1; % HARDCODED
J_lowpass=J.*low_pass;
j_lowpass = ifft2(fftshift(J_lowpass));
figure;
subplot(1,3,1);
imshow(cropped_array)
title('Original Image')  
subplot(1,3,2);
imshow(abs(j_lowpass),[])
title('Flat Field + Low Pass')  
subplot(1,3,3);
imshow(log(abs(J_lowpass)),[]);
title('Fourier Domain')  

B = rescale(abs(j_lowpass));
%% non local filter
imshow(B)
roi = drawrectangle;
roi_position = roi.Position;
patch = imcrop(B,roi_position);
patchSq = patch.^2;
edist = sqrt(sum(patchSq,3));
patchSigma = sqrt(var(edist(:)));
DoS = 1.5*patchSigma;
%%
filtered_im = imnlmfilt(B,'DegreeOfSmoothing',DoS);
imshow(filtered_im)
hold on
rectangle('Position', roi_position )

%% binarize
B_f=1-filtered_im;
BW=imbinarize(B_f,"adaptive",'ForegroundPolarity','bright','Sensitivity',0.58); % HARDCODED
%BW = imfill(BW,"holes");

% BW=imbinarize(B_f,0.47);

figure;
subplot(1,2,1);
imshow(BW,[])
se = strel('rectangle',[7 12]);
sk = imclose(BW,se);
se = strel('rectangle',[12 7]);
sk = imclose(sk,se);
se = strel('line',3,45);
sk = imclose(sk,se);
sk = imopen(sk,se);
se = strel('line',3,-45);
sk = imclose(sk,se);
sk = imopen(sk,se);
subplot(1,2,2);
imshow(sk,[])

%% Skeleton
out = bwskel(sk, 'MinBranchLength',15); % HARDCODED
figure;
subplot(1,2,1);
imshow(input,[])
subplot(1,2,2);
imshow(out)
%% overlay bead
[columnsInImage, rowsInImage] = meshgrid(1:size(B,1), 1:size(B,1)); 
circlePixels = (rowsInImage - centers(2)).^2 + (columnsInImage - centers(1)).^2 <= radii(1).^2;
inner_circle=(rowsInImage - centers(2)).^2 + (columnsInImage - centers(1)).^2 <= (radii(1)-1).^2;
skel=logical(out+circlePixels)-inner_circle;
figure;
imshow(labeloverlay(B,skel,'Transparency',0,'Colormap','autumn'))

% https://matlab.fandom.com/wiki/FAQ#How_do_I_create_a_circle.3F
%% labeling
close all;
[L,obg_num]=bwlabel(skel,8);
b=L(floor(centers(2)),floor(centers(1)+radii(1)));
L=(L==b);
figure;
subplot(1,2,1);
imshow(1-L,[])
title('Skeleton')
subplot(1,2,2);
imshow(labeloverlay(B,L,'Transparency',0,'Colormap','autumn'))
title('Overlay')
%% measure length
sprouts=L-circlePixels;
sprouts=sprouts>0;

[sprouts_L,num_sprouts]=bwlabel(sprouts,8);
imshow(1-sprouts_L,[])

lengths=zeros(1,num_sprouts);

for i=1:num_sprouts
    lengths(i)=sum(sprouts_L==i,'all');
end
lengths=lengths./scale;
average_length=mean(lengths);
total_length=sum(lengths,"all");

%% Export results
image_num=regexp(myDir,'\','split');
image_num=image_num(end);
filename=convertCharsToStrings(myFiles(slice_selected).name);
filename=regexp(filename,'_','split');
date=filename(1);
experiment=filename(2);
slice=slice_selected-1;

T = table(date,experiment,image_num,slice,diameter,num_sprouts,lengths,average_length,total_length,sensitivity_bin);
writetable(T,'testing_pipeline_2.xlsx','WriteMode','Append');

output_image_file = convertCharsToStrings(myFiles(slice_selected).name);
[filepath,name,ext] = fileparts(output_image_file);
output_image_file = 'Skeleton_'+ date + '_'+ experiment + '_'+ image_num + '_'+ 'z' + slice + '.png';
imwrite(logical(L),output_image_file)
output_mask_file = 'Overlay_'+ date + '_'+ experiment + '_'+ image_num + '_'+ 'z' + slice + '.png';
imwrite(imfuse(L,B,'blend','Scaling','joint'), output_mask_file);