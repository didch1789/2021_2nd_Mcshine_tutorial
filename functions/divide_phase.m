function [out, crit] = divide_phase(in, varargin)

drawnow = false;

for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'drawnow'
                drawnow = true;
        end
    end
end

PC1t = in(:, 1);
dPC1t = diff([PC1t(1);PC1t]);

lowcrit = prctile(PC1t, 25);
highcrit = prctile(PC1t, 75);

low = PC1t <= lowcrit;
rise = (PC1t >= lowcrit) & (PC1t <= highcrit) & dPC1t > 0;
high = PC1t >= highcrit;
fall = (PC1t >= lowcrit) & (PC1t <= highcrit) & dPC1t < 0;

if drawnow
    scatter(find(low), PC1t(low), 'filled');hold on;scatter(find(rise), PC1t(rise), 'filled');
    scatter(find(high), PC1t(high), 'filled');scatter(find(fall), PC1t(fall), 'filled');
end

out.low = low;
out.rise = rise;
out.high = high;
out.fall = fall;
crit.low = lowcrit;
crit.high = highcrit;

end


