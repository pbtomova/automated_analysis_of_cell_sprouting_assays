function [centers,radii]  = detect_beads(x)
    
    J = imflatfield(x,60); % HARDCODED
    log_im=(3*log(1+im2double(J)));
    
    figure;
    subplot(1,2,1);
    imshow(x);
    title('Original image')
    
%     subplot(2,2,3);
%     imhist(x,100) 
%     subplot(2,2,4);
%     imhist(log_im,100)
    filtered_im=filter2(fspecial('average',3),log_im);
    
    subplot(1,2,2);
    imshow(filtered_im);
    title('Flatfield Correction + Log Transformation image')

    BW=imbinarize(filtered_im,"adaptive",'ForegroundPolarity','dark');
    
    se = strel('disk',5,0);
    BW = imopen(BW,se);
    BW = imopen(BW,se);


%     se = strel('line',10,0);
%     BW = imopen(BW,se);
%     BW = imopen(BW,se);
% 
%     se = strel('line',10,90);
%     BW = imopen(BW,se);
%     BW = imopen(BW,se);
    
     % figure

    [centers,radii] = imfindcircles(BW,[30 100],'ObjectPolarity','dark', ...
        'Sensitivity',0.91);
    figure;
    imshow(log_im)
    h = viscircles(centers,radii);

    
end