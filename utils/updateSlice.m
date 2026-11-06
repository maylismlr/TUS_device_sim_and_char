function updateSlice(src, hImg, pressure)
% Little function to add a slicer for visualisations
%
%   author: Maylis Mueller

    idx = round(src.Value);
    hImg.CData = pressure(:,:,idx);
    title(['Slice: ', num2str(idx)]);
end