function dice = get_DSC(mask_A, mask_B)
% get dice similarity coefficient between two masks
%
%   author: Maylis Mueller

    denominator = nnz(mask_A) + nnz(mask_B);

    if denominator > 0
        dice = ...
            2 * nnz(mask_A & mask_B) / denominator;
    else
        dice = NaN;
    end
end
