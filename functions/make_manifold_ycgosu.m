function manifold = make_manifold_ycgosu(phasemat, onetPC, refTime)

if nargin < 3
    refTime = 30;
end

for ph = 1:size(phasemat, 2)
    phidx = extractones(phasemat(:, ph));
    meanTRs = round(mean(cellfun(@(x) size(x, 2), phidx)));
    interped = [];  iidx = 1;
    for phi = 1:numel(phidx)
        if (numel(phidx{phi}) < 2)
            continue
        else
            interped(iidx, :, :) = interp1(phidx{phi}, onetPC(phidx{phi}, :), linspace(phidx{phi}(1), phidx{phi}(end), refTime));
            iidx = iidx + 1;
        end
    end
    manifold(ph, :, :) = squeeze(mean(interped, 1));
end
    

end