function M = focal_metrics(P, kgrid)
% Peak pressure and sub-grid focal position.
%
% The pressure maximum is searched over the complete computational domain.
% Parabolic sub-grid refinement is performed only if neighbouring samples
% exist on both sides.
%
%   author: Maylis Mueller

    [M.pmax, li] = max(P(:));
    [i, j, k] = ind2sub(size(P), li);
    M.idx = [i, j, k];
    
    dxs = 0;
    dys = 0;
    dzs = 0;
    
    if i > 1 && i < size(P,1)
        dxs = parab(P(i-1,j,k), P(i,j,k), P(i+1,j,k));
    end
    
    if j > 1 && j < size(P,2)
        dys = parab(P(i,j-1,k), P(i,j,k), P(i,j+1,k));
    end
    
    if k > 1 && k < size(P,3)
        dzs = parab(P(i,j,k-1), P(i,j,k), P(i,j,k+1));
    end

     M.pos = [kgrid.x_vec(i) + dxs*kgrid.dx, ...
             kgrid.y_vec(j) + dys*kgrid.dy, ...
             kgrid.z_vec(k) + dzs*kgrid.dz];
end


function d = parab(y1, y2, y3)
% Three-point parabolic interpolation around a discrete maximum.

    den = y1 - 2*y2 + y3;

    if den == 0

        d = 0;
        return;
    end

    d = 0.5 * ...
        (y1 - y3) / den;
end