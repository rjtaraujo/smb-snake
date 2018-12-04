function H = handRegionGrowing(D)
% This function employs region growing to obtain a preliminary segmentation
% of the hand region in a depth image.
%   inputs: I - the depth image of the hand
%           
%   outputs: H - a preliminary mask of the hand region

        sorted_values = unique(sort(D(:)));
        depth_corr = D;
        depth_corr(depth_corr > sorted_values(end-2)) = sorted_values(end-2);
        depth_corr(depth_corr < sorted_values(3)) = sorted_values(3);
        depth_corr = (depth_corr - min(depth_corr(:))) / (max(depth_corr(:)) - min(depth_corr(:)));
        
        % Entropy of depth information
        entr = entropyfilt(depth_corr,ones(5,5));
        entr = (entr - min(entr(:))) / (max(entr(:)) - min(entr(:)));
    
        % Cleaning the depth data
        comb = (1 - depth_corr) .* (1 - entr);
    
        % Finding a seed inside the hand
        K = fspecial('average', 25);
        avg_I = imfilter(comb, K);
        [~,ind] = max(avg_I(:));
        ind = ind(1);
        [row, col] = ind2sub(size(avg_I), ind);

        this_RG = regiongrowing(comb, row, col, 0.23);        
        H = im2double(imclose(this_RG, strel('disk',2)));

end

