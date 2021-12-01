function coordinate = mni2cor(mni, T)

if isempty(mni)
    coordinate = [];
    return;
end

if nargin == 1
	T = ...
        [-2     0     0    92;...
         0     2     0  -128;...
         0     0     2   -74;...
         0     0     0     1];
end

coordinate = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]*(inv(T))';
coordinate(:,4) = [];
coordinate = round(coordinate);

end


