%% Mcshine 2019.
% Prerequisite. 
%   
% Differences in the following code. 
%
%   1. Used resting state (same TR (720ms) but different total length (1200TRs)
%   2. Preprocessing was done in a bit different way.
%       - Mcshine used minimally preprocessed data (upto field & mo. co.)
%         and did the following preprocessing (similar with ours, but only low-pass)
%       - Following haved used ICA-FIX + Band Pass filtered (.008 ~ .1Hz)
%   3. Same processes.
%       - No smoothing.
%       - Used same atlas. (Is it right to use Gordon atlas? It's surface atlas...)
basedir = '/Users/jungwookim/Dropbox/github/2021_2nd_classes/2021_Woo/Mcshine_tutorial';
addpath(genpath(fullfile(basedir, 'functions')))
cd(basedir)

%% 1. Calculate spatial PCs.
load(fullfile(basedir, 'roi_extracted.mat')); % roiVals
% Loading concatenated subjects' resting state data.
% Each cell contains subject data(20 subjs) / 1200 TRs (1200 * 0.72 : 14min 30sec)

[coeffs, ~, ~, ~, expls] = pca(cell2mat(roiVals), 'NumComponents', 5);
bar(expls(1:5)) % First PC accounts for 30% while top 5 explains 46%. / orig paper: 68%

for i = 1:numel(roiVals)
    subcoeff{i} = pca(roiVals{i}, 'NumComponents', 5);
end

% Are they similar acorss individuals?
for i = 1:numel(subcoeff)
    subPCwholPC(:, :, i) = corr(coeffs, subcoeff{i});
end
subplot 211;
imagesc(mean(subPCwholPC, 3));colorbar;
title({'PC cos similarity averaged across Subjs'})
subPC1s = cell2mat(cellfun(@(x) x(:, 1), subcoeff, 'UniformOutput', false));
subplot 212;
imagesc(corr(subPC1s));colorbar;
title('PC1 similarity across Subjs')
% seems so!
set(gcf, 'Position', [0 0 470 580])

%% 1-1. Visualize Spatial PCs. (can skip)
% Require canlab, cocoan, & yangcho_gosu toolbox.
% Unwrap & run following three 3 command lines.
% system('git clone https://github.com/canlab/CanlabCore')
% system('git clone https://github.com/cocoanlab/cocoanCORE') 

% Saving in nifti files. (Already in ...
% fullfile(basedir, 'results', 'SpatialPCs')
% outorig = map2original(coeffs);
% for i = 1:size(outorig.dat, 2)
%     temp = outorig;
%     sPC = outorig.dat(:, i);
%     temp.dat = sPC .* (sPC > prctile(sPC, 95) | sPC < prctile(sPC, 5));
%     temp.fullpath = fullfile(basedir, 'results', 'SpatialPCs', sprintf('spatialPC%d_thr5per.nii', i));
%     write(temp, 'overwrite');
%     no info. of how they have thresholded the results...!
%     so arbitrarily used 5 % threshold.
% end

% Plotting... (Already in "
% fullfile(basedir, 'results', 'SpatialPCs') 
close all;
for i = 1:5
    figure(i)
    brain_activations_display(region(which(sprintf('spatialPC%d_thr5per.nii', i))), 'depth', 2, 'surface_only', ...
        'surface_all')
%     saveas(gcf, fullfile(basedir, 'results', 'SpatialPCs', sprintf('spatialPC%d_thr5per', i)), 'png');
%     close all
end

% close all
% numPC = 3; % choose between 1~5
% radialplot_network_overlap(which(sprintf('spatialPC%d_thr5per.nii', numPC)));
% Seems bit different...!

%% 2. Calculating tPCs 
cd(fullfile(basedir))
clearvars -except roiVals coeffs basedir;clc;

for i = 1:numel(roiVals)
    onesubtPC = roiVals{i} * coeffs;
    tPCs{i} = onesubtPC;
end
cattPC1s = cat(3, tPCs{:});
sccattPC1s = squeeze(cattPC1s(:, 1, :));
imagesc(corr(sccattPC1s));colorbar % substantially low between subject correlation...

% Matching process could be done!!
% Since we're assuming that all indiviudals have same latent, oscillaotory
% fluctuations, we focus on that by using *QPP and *CAP...
% But we're interested in MANIFOLD structure! so we're going to 
% draw MANIFOLD structure in each participant, and try to 
% see if similar pattern occur..
% *QPP: Quasi-Periodic Pattern
% *CAP: Co-activation Pattern analysis.

plot(tPCs{1}(:, 1)); % plot of one participant tPC1. still quite noisy though...!

%% Additional processing of tPC. (needs discussion)
% only extracting frequency with high power!
% Original paper normalized the data within task,
% to match the amplitude of the whole time series.
% Since we have no such time points, I decided to
% use only top 99 percent amplitude frequencies to
% get the "peak normalization effect" so that we can
% extract proper time points of each phase...

sF = 1/0.72; % sampling frequency
T = 0.72; % TR 
L = 1200; % Total length
% t = (0:L-1)*T; % Total time in sec.
% f = sF*(0:L/2)/L; % for freq domain.

for i = 1:numel(tPCs)
    onetPC = tPCs{i};
    for nP = 1:size(onetPC, 2)
        fftPC = fft(onetPC(:, nP));
        P2 = abs(fftPC/L);
        P1 = P2(1:L/2);
        realYtPC = real(fftPC);  imagYtPC = imag(fftPC);
        indx = abs(fftPC) > prctile(abs(fftPC), 99); % Only using top 1% power of the spectrum.
        recontPC{i}(:, nP) = ifft((realYtPC .* indx) + (imagYtPC .* indx) * 1i);
    end 
end
subplot 211
divide_phase(tPCs{1}(:, 1), 'drawnow');hold on
title('Original tPC1')
subplot 212
divide_phase(recontPC{1}(:, 1), 'drawnow');
title('Reconstructed tPC1')
set(gcf, 'Position', [0 0 470 580])


%% Reconstructing tPCs and making manifolds...
% This part is bit suspiscious...
% In order to make a manifold, they had to make same number of length for
% each trajectory. (i.e., one cycle has 4 trajectory. and one trajectory (e.g. low) 
% has different length in the next cycle.) So they just used interpolation
% to match the length between cycles. and then average them across all
% cycles. 
manifolds = {};

for i = 1:numel(tPCs)
    onetPC = tPCs{i};
    phaseout = divide_phase(recontPC{i}(:, 1)); 
    phasemat = [phaseout.low phaseout.rise phaseout.high phaseout.fall];
    manifolds{i} = make_manifold_ycgosu(phasemat, onetPC);
end

%% Group manifold.
cmaps = [0 0 195;193 0 8;252 59 2;74 166 250]./255;
avgsubmanifold = mean(cat(4, manifolds{:}), 4);

fig1 = figure(1);
for j = 1:4
    scatter3(avgsubmanifold(j, :, 1), avgsubmanifold(j, :, 2), avgsubmanifold(j, :, 3), ...
        40, 'filled', 'MarkerFaceColor',cmaps(j, :)); hold on;
    xlabel('tPC1');ylabel('tPC2');zlabel('tPC3')
end  
% bump attractor?

fig2 = figure(2);
for i = 1:numel(manifolds)
    submanifold = manifolds{i};
    for j = 1:4
        scatter3(submanifold(j, :, 1), submanifold(j, :, 2), submanifold(j, :, 3), ...
            40, 'filled', 'MarkerFaceColor',cmaps(j, :)); hold on;
        xlabel('tPC1');ylabel('tPC2');zlabel('tPC3')
    end  
end

%% Loading on Allen brain atlas.
ColorCodeAccordingToExcitInHibit = true;
if ColorCodeAccordingToExcitInHibit
    cmaps = [.8 .2 .2;.2 .2 .8;.2 .2 .8;.8 .2 .2;.8 .2 .2;.2 .2 .8;.8 .2 .2];
else
    cmaps = ...
        [0.8941    0.1020    0.1098
        0.2157    0.4941    0.7216
        0.3020    0.6863    0.2902
        0.5961    0.3059    0.6392
        1.0000    0.4980         0
        1.0000    1.0000    0.2000
        0.6510    0.3373    0.1569];
end
load roi_extracted.mat
load NT_timeseries_smth.mat

% Processing procedures of Allen Brain atlas was unclear...
% Allen Brain Atlas:
%   samples: Samples within single specimen. 
%   probes: contains which gene expression was measured.
%   There are several probes per sample.
%   (sample can be regarded as location and probes indicate
%    information (of gene expression) within region.)
%   Only around 1000 voxel information was provided...
%   Since voxel level map might be inappropriate for our usage,
%   I smoothed with 4mm kernel and used it as receptor map.
%   Resting state * Receptor map = one time series.
%   Compare the resulting time series and tPCs...
%       Seems direct comparison between spatial PCs and receptor map
%       is more ideal... 

for i = 1:numel(tPCs)
    subNTs = squeeze(NTtimes(i, :, :));
    subtPC = recontPC{i};
    
    for nt = 1:size(subNTs, 2)
        subcorrsPC1(i, nt) = corr(subtPC(:, 1), subNTs(:, nt));
        subcorrsPC2(i, nt) = corr(subtPC(:, 2), subNTs(:, nt));
        subcorrsPC3(i, nt) = corr(subtPC(:, 3), subNTs(:, nt));
    end
end

for i = 1:19
    subplot(5, 4, i)
    for j = 1:7
        scatter(subcorrsPC1(i, j)', subcorrsPC2(i, j)',...
            'filled', 'MarkerFaceColor', cmaps(j, :));hold on
    end
end
legend(Geneinfos, 'Interpreter', 'none', 'Location', 'southeastoutside')
% No consistent results found...!


























