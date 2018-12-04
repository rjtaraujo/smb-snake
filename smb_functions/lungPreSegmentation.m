function L = lungPreSegmentation(I)
% This function employs Otsu's threshold to obtain a preliminary segmentation
% of the lung region in a CT scan.
%   inputs: I - the CT scan image
%           
%   outputs: L - a preliminary mask of the lungs region

    th_body = graythresh(I(I>0));
    temp = logical(im2bw(I,th_body));
        
    CC = bwconncomp(temp,4);
    stats = regionprops(CC,'PixelIdxList','Area');
    areas = [stats.Area];
    
    idx = find(areas==max(areas));
    body = zeros(size(I));
    body(stats(idx).PixelIdxList) = 1;
    
    L = im2double(bwareaopen(imfill(body,4,'holes') - body, 2500, 4));        
end

