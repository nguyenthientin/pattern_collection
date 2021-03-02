"""
 @author Tin T. Nguyen
 @email t.t.nguyen-3@tudelft.nl
 @create date 2020-11-03
 @modify date 2020-11-03
 @desc [description]
""" 

import numpy as np
import warnings
from constants import SPEED_NaN, FLOW_NaN

def segment_congestion_pattern(speed, raw_speed, flow=None, raw_flow=None, 
                                congestion_naive_mask=None,
                                ):
    """SEGMENT_CONGESTION_PATTERN segment a congestion patterns to identify wide moving jam,
    homogeneous region, bottlenecks.

    Args:
        speed ([type]): [description]
        raw_speed ([type]): [description]
        flow ([type], optional): [description]. Defaults to None.
        raw_flow ([type], optional): [description]. Defaults to None.
        congestion_naive_mask ([type], optional): [description]. Defaults to None.
    """
    # handle NaN values in speeds
    if np.any(np.isnan(speed)):
        warnings.warn('NaN(s) found in speed. Replace those with maximum speed.'
                       + 'The result might be unreliable')
    np.max()
    nanMask = speed == speedNaN;
    speed(nanMask) = -1;
    speed(nanMask) = max(speed(:));
    pData.speed = flipud(speed);
%% false detectors
if ~isempty(pData)
    nanMaskRaw = pData.raw.speed == 99999;% isnan(pData.raw.speed);
    s = mean(nanMaskRaw, 2);
    falsedetind = find(s >= 0.7);
    pData.raw.x(falsedetind) = [];
    pData.raw.speed(falsedetind, :) = [];
    pData.raw.flow(falsedetind, :) = [];   
end
%%
%speed, flow, centercongmask, freeMask, outdir, fn)
tic;
% figure('units','normalized','outerposition',[0 0 1 1]) % for all plots
% figure
% axF = axes;
% hold on
%%
% make sure speeds are double numbers
speed = double(speed);
%%
EXTREME_CONGESTED_SPEED = 20;
EXTREME_CONGESTED_SPEED = 10;
FREE_SPEED = 90;
spdROI = speed;
%% smooth
spdBlur = imgaussfilt(spdROI, 0.5);
%%
congMask = cpc_hp_estimatecongestedmask(speed, 'contour_builtin', 0);
if ~any(congMask(:))
    info.r = [];
    return
end
% focus on the center patterns
congcenter = flipud(pData.congmask);
[cl, ncl] = bwlabel(congMask);
for i = 1:ncl
    ap = cl == i;
    cin = sum(congcenter(:) & ap(:));
    if cin < sum(ap(:)) / 3
        congMask(ap) = 0;
    end
end
% c = cl(pData.congmask);
% cunik = unique(c); cunik(cunik == 0) = [];
% congMask = ismember(cl, cunik);

% ut = speed ./ (max(speed(:)));
% congMask = ~imbinarize(ut);
debug = 0;
if debug
    figure; imagesc(congMask);
    figure; imagesc(speed);
end
freeMask = ~congMask;
% shrunkenFreeMask = imerode(freeMask, ones(5));
shrunkenFreeMask = freeMask;
%% edge & gradient
% [spdEdge] = edge(spdROI, 'canny', 0);
[spdEdge] = edge(spdROI, 'canny', [], 0.7);
% e = edge(spdROI, 'canny');
% spdEdge = e | spdEdge;
% [spdEdge] = edge(spdROI, 'Sobel', 0);
% gradient
[Gmag, Gdir] = imgradient(spdBlur);
% Gmag = pr_domultiscalegradient(spdBlur, 2);
% watershed
% Gmag = imgaussfilt(Gmag, 1);
%% clean & normalize Gmag
Gmag = Gmag / max(Gmag(:));
% Gmag = imclose(Gmag, ones(1));
[gRef, g, fi, rawcongmask, btnInfo, asmbtnline] = btn_detectbtnfromrawdata(pData);
% test11raw;
% test11_refcongmask;
doOld = 1;
if doOld
    Gmag(spdEdge == 1) = 1;
    Gmag = imhmin(Gmag, 0.015);
    spdCong = spdROI < EXTREME_CONGESTED_SPEED;
    Gmag(spdCong) = min(Gmag(spdCong));
    % spdFree = spdROI >= 100; %FREE_SPEED;
    % spdFree = imerode(spdFree, ones(5));
    Gmag(spdEdge == 1) = 1;
    % Gmag(spdFree) = min(Gmag(:));
    shrunkenFreeMask = imerode(~congMask, ones(10));
    Gmag(shrunkenFreeMask) = min(Gmag(:));
    % Gmag(Gmag < 0.03) = 0;
    %% watershed
    % WGmag = watershed(Gmag, 4);
    WGmag = watershed(Gmag);
    %% remove free traffic regions
    % freeMask = data_hp_getpatternmask(speed, freeMask);
    % freeMask = ~spdFree; % this is temporary
    % WGmagFreeIds = unique(WGmag(~freeMask));
    WGmagFreeIds = unique(WGmag(shrunkenFreeMask));
    WGmag(ismember(WGmag, WGmagFreeIds)) = 0;
    % WGmag = double(WGmag) .* freeMask;
    WGmag = pr_cleansegmentation(WGmag, congMask);
    %% retain only congested regions

    %% post processing
    % WGmag = pr_postprocess(WGmag, spdROI, 80, 0.5);
    info.WId = WGmag;
    %% represent basins as mean speed
    spdBasins = pr_calspeedbasins(WGmag, speed);
    info.WSpd = spdBasins;
    %% merging
    % adapt speed and watershed basins
    shed = spdBasins; shed(spdBasins ~= 120) = 0;
    edgeNotShedInd = spdEdge & (~shed);
    spdEdge(edgeNotShedInd) = 0;
    WGmagMerge = pr_mergeregions(WGmag, spdROI, 120, spdEdge);
    %% represent merging basins by mean speed
    spdMergeBasins = pr_calspeedbasins(WGmagMerge, speed);
    %% remove shed lines
    [WGmagMerge, spdMergeBasins] = cpc_removeshedlines(...
                                        WGmagMerge, spdMergeBasins, speed);
    info.WMId = WGmagMerge;
    info.WMSpd = spdMergeBasins;
    %% save binarization
    % subplot(1,2,1); imagesc(spdMergeBasins);
    % subplot(1,2,2); imagesc(spdBlur);
    % saveas(gcf, [pSavePathName, '.jpg']);
    % return
    %% find WMJ candidates
    if ~isempty(WGmagMerge)
        [wmjSpd, wmjIds] = ...
            cpc_extractdisturbances(spdBlur, WGmagMerge, spdMergeBasins);
            %cpc_extractwmj(spdBlur, WGmagMerge, spdMergeBasins);
    end
    info.wmjId = wmjIds;
    wmjIdList = unique(wmjIds(:));
    wmjIdList(wmjIdList == 0) = [];
    %% 
    btnconglabel = btn_raw_separatecongestion(congMask, asmbtnline);
    nj = max(btnconglabel(:));
    btnwmjlabel = btnconglabel;
    btnwmjlabel(wmjIds ~= 0) = nj + 1;
    %figure; imagesc(btnwmjlabel)
    %% create a (graph-based) representation of this pattern
    info.r = btn_gengraphrepresentation2(congMask, btnconglabel, wmjIds);
    return
    %% 
    wmIds = WGmagMerge;
    wmSpd = spdMergeBasins;
    spdIdMap = cpc_hp_generateidvaluemap(wmIds, wmSpd);

    congMask = imclose(wmIds ~= 0, ones(3)); % mask of congestion
    % congMask = wmIds ~= 0;
    connComps = bwconncomp(congMask); % connected components

    MIN_SIZE = 50;
    btnIds = zeros(size(speed));
    btnLocs = [];
    btnVis = speed;
    % stationary bottleneck only
    staIds = zeros(size(speed));
    % staLocs = [];
    btnInfo = [];
    for iC = 1:connComps.NumObjects
        compInd = connComps.PixelIdxList{iC};
        if length(compInd) < MIN_SIZE
            continue % ignore if the partition is too small
        end
        isCenterCong = centercongmask(compInd);
        if ~any(isCenterCong)
            continue
        end
        % get the corresponding congestion partition
        pIds = zeros(size(speed));
        pIds(compInd) = wmIds(compInd);
        pSpd = 120 * ones(size(speed));
        pSpd(compInd) = wmSpd(compInd);


    %     % remove wmj 
    %     rInd = ~ismember(pIds, wmjIdList); % indices of the remain
    %     rSpd = pSpd;
    %     rSpd(~rInd) = 120;
    %     rIds = pIds;
    %     rIds(~rInd) = 0;
    %     % remove intermediate regions
    % 
    %     [rIds, rSpd] = cpc_filterintermediateregion(...
    %                         rIds, rSpd, wmjIds);

        % bottleneck
        HOMOGENEITY_DT_MIN = 40; % 20 mins
        homo = pSpd < 30;
        homo = imopen(homo, ones(1, HOMOGENEITY_DT_MIN));
        metadata = struct('speed', speed, 'spdBlur', spdBlur, ...
                          'homo', homo);
        [rgsId, ~, ~, locs, isTimeout, rgsBtnVis, rgsBtnInfo] = btn_detectbottlenecks(...
                                pIds, pSpd, wmjIds, wmjSpd, pIds, metadata);
    %     [rgsId, ~, ~, locs] = btn_detectbottlenecks(rIds, rSpd, wmIds);
        if isTimeout
            disp('timeout')
            disp(outdir)
        end
        if isempty(rgsId)
            continue
        end

        btnLocs = [btnLocs; locs];
        btnInfo = [btnInfo rgsBtnInfo];
        % add to the static btn
        sId = rgsId;
        sId(sId ~= 0) = sId(sId ~= 0) + max(staIds(:));
        staIds = staIds + sId;

        % add back wmjs
        pMask = pIds ~= 0;
        rMask = rgsId ~= 0;
        jMask = pMask & ~rMask;
        % rgsId = btn_addBackWMJ(rgsId, jMask);

        % add to the final bottleneck
        pBtnId = rgsId;
        pBtnId(pBtnId ~= 0) = pBtnId(pBtnId ~= 0) + max(btnIds(:));
        btnIds = btnIds + pBtnId;

        % visualization
        btnVis(logical(rgsBtnVis)) = 0;
    end
    %% debug
    figure; 
    subplot(1,3,1); imagesc(wmSpd)
    subplot(1,3,2); imagesc(btnIds);
    subplot(1,3,3); imagesc(btnVis);
end
%% plot
figure;
subplot(121); imagesc(speed); colormap('gray')
uvis = speed;
uvis(g) = 0;
subplot(122); imagesc(uvis)
% subplot(122); imagesc(uvis); colormap(flipud(jet(max(speed(:)))))
%% save
ISSAVE = p.Results.save;
if ~ISSAVE
    return
end
% save figure
outdir = p.Results.outdir;
fn = p.Results.filename;
outdirimg = fullfile(outdir, 'img');
if ~exist(outdirimg, 'dir')
    mkdir(outdirimg)
end
saveas(gcf, fullfile(outdirimg, [fn, '.jpg']));
% save btninfo
outdirinfo = fullfile(outdir, 'info');
if ~exist(outdirinfo, 'dir')
    mkdir(outdirinfo)
end
btnLocs = find(any(gRef, 2));
btnInfo = [];
for i = 1:length(btnLocs)
    btnInfo(i).id = i;
    btnInfo(i).loc = btnLocs(i);
end
save(fullfile(outdirinfo, [fn, '_btnInfo.mat']), 'btnInfo')
close all
% figlist = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(figlist)
%     fighandle = figlist(iFig);
%     figname = [pSavePathName, '_part', num2str(iFig), '.jpg'];
%     figure(fighandle.Number)
%     saveas(gcf, figname);
%     close
% end
% saveas(gcf, [pSavePathName, '.jpg']);
% print(pSavePathName, '-dpng', '-r0');
% close all
return
%%
info.btnId = btnIds;
info.btnLocInd = btnLocs;
% %% extract WMJ candidates
% if ~isempty(WGmagMerge)
%     [wmjSpdMergeBsins, wmjSpdIds] = ...
%         cpc_extractwmj(spdBlur, WGmagMerge, spdMergeBasins);
% end
% %% save
% % bwRgs = wmjSpdIds ~= 0;
% % imwrite(bwRgs, [pSavePathName, '.jpg']);
% % 
% % WMJCand = bwRgs;
% % Congestion = WGmagMerge ~= 0;
% % 
% % save([pSavePathName, '.mat'], 'WMJCand', 'Congestion')
% 
% %% get the remainning
% wmjIds = unique(wmjSpdMergeBsins(:));
% wmjIds(wmjIds == 120) = [];
% remainPatternInd = ~ismember(spdMergeBasins, wmjIds);
% remainPattern = spdMergeBasins;
% remainPattern(~remainPatternInd) = 120;
% remainPatternId = WGmagMerge;
% remainPatternId(~remainPatternInd) = 0;
% % remove intermediate regions
% 
% [remainPatternId, remainPattern] = ...
%     cpc_filterintermediateregion(remainPatternId, remainPattern, wmjSpdIds);
% %% not removing wmjs
% % remainPattern = spdMergeBasins;
% % remainPatternId = WGmagMerge;
% %% bottleneck
% spd2Value = remainPattern;
% spd2Id = remainPatternId;
% [rgsId, rgsRange, elapsedTime] = btn_detectbottlenecks(spd2Id, spd2Value);
%% visualize
spdRange = [0 150];
pathParts = strsplit(pSavePathName, filesep);
% figure('units','normalized','outerposition',[0 0 1 1])
nRP = 2; % plot: number of rows
nCP = 3; % plot: number of columns
warning('off','all')
subplot(nRP,nCP,1); imagesc(spdROI); title(['speed - ', pathParts{end}]);
subplot(nRP,nCP,4); imagesc(flow);
subplot(nRP,nCP,2); imagesc(wmjSpd); title('Wide moving jams');
subplot(nRP,nCP,5); imagesc(spdMergeBasins); title('Merge - Speed'); truesize
%
btnMask = btnIds ~= 0;
removedCongMask = congMask & ~btnMask;
btnIds(removedCongMask) = max(btnIds(:)) + 1; % to display the rest congestion in a different color
subplot(nRP,nCP,3); imagesc(btnIds); %imagesc(rgsId);
% subplot(2,3,6); imagesc(spdEdge);
subplot(nRP,nCP,6);
spdWBtnLocs = speed;
spdWBtnLocs(btnLocs, :) = 0;
imagesc(spdWBtnLocs)
hold on
% boundary
nBtn = numel(btnLocs);
colMap = jet(nBtn);
btnBoundMethod = 'wch';
for i = 1:numel(btnLocs)
    [btnConvBound, boundMask] = ...
        btn_getboundary(btnLocs(i), i, btnBoundMethod, congMask, staIds, wmjIds);
               % btn_getboundary_v2(btnLocs(i), i, btnIds, staIds, congMask);
%     spdWBtnLocs(boundMask) = 0;
    plot(btnConvBound(:,2), btnConvBound(:,1), ...
        'Color', colMap(i,:), ...
        'LineWidth', 2)
end
hold off
% imagesc(spdWBtnLocs)
%% save
ISSAVE = false;
if ~ISSAVE
    return
end
saveas(gcf, [pSavePathName, '.jpg']);
% print(pSavePathName, '-dpng', '-r0');
close all
%%
time = toc;
fprintf('%2.2f seconds\n', time);
end
function WGmagPost = pr_postprocess(WGmag, speed, thr, pTile)
WGmagPost = WGmag;
nBasins = max(WGmag(:));
% FREE_BASIN_SPEED = 65;
%% filt out un-cong regions
for i = 1:nBasins
    bInd = WGmag == i;
    bSpeed = speed(bInd);
    bSize = sum(bInd(:));
    unCongArea = sum(bSpeed > thr);
    if unCongArea / bSize > pTile
        WGmagPost(bInd) = 0;
    end
end
end
%%
function WGmagPost = pr_filterspeed(WGmag, speed, thr)
WGmagPost = WGmag;
nBasins = max(WGmag(:));
% FREE_BASIN_SPEED = 65;
%% filt out above-threshold regions
for i = 1:nBasins
    bInd = WGmag == i;
    bSpeed = speed(bInd);
    if isempty(bSpeed)
        continue
    end
    if bSpeed(1) > thr
        WGmagPost(bInd) = 0;
    end
end
end
%%
function spdBasins = pr_calspeedbasins(WGmag, speed)
spdBasins = speed;
bsIds = unique(WGmag((WGmag > 0)));
nBs = numel(bsIds);
for iBs = 1:nBs
    bsId = bsIds(iBs);
    bsInd = WGmag == bsId;
    bsSpd = speed(bsInd);
    spdBasins(bsInd) = mean(bsSpd);
end
spdBasins(WGmag == 0) = 120;
end
function multigrad = pr_domultiscalegradient(image, n)
multigrad = zeros(size(image));
for i = 1:n
%     b = strel('disk', 2*i+1);
%     bPre = strel('disk', 2*i-1);
    b = ones(2*i+1);
    bPre = ones(2*i-1);
    grad = imerode((imdilate(image, b) - imerode(image,b)), bPre);
    multigrad = multigrad + grad;
end
multigrad = multigrad / n;
end
%%
function finalRgsSeg = pr_mergeregions(rgsSeg, rgsValue, thresh, edge)
% shedlines
SHEDLINE_ID = 0;
shed = rgsSeg == SHEDLINE_ID;
% get interesting regions
rgsROI = rgsValue < thresh;
rgs.Seg = rgsSeg;
rgs.SegValue = rgsValue;
% rgs.Ind = rgsSeg < thresh;
rgs.Ids = unique(rgsSeg(rgsROI));
rgs.Ids(rgs.Ids == SHEDLINE_ID) = []; % remove shed lines id
rgs.nIds = numel(rgs.Ids);
% order this Ids in increasing speed
spds = zeros(1, rgs.nIds);
for iId = 1:rgs.nIds
    id = rgs.Ids(iId);
    regionSpd = rgsValue(rgsSeg == id);
    spds(iId) = mean(regionSpd);
end
[~, spdOrder] = sort(spds, 'ascend');
rgs.Ids = rgs.Ids(spdOrder);
%
rgs.isProcessed = zeros(rgs.nIds, 1);
for iId = 1:rgs.nIds
    if rgs.isProcessed(iId)
        continue
    end
    % grow this region (if possible)
%     edge = imdilate(edge, ones(2));
%     edge = edge & shed;
%     disp(rgs.Ids(iId));
    rgs = pr_growandmergeregion(rgs, iId, edge);
    rgs.isProcessed(iId) = 1;
%     imagesc(rgs.Seg);
end
finalRgsSeg = rgs.Seg;
end
%%
function rgs = pr_growandmergeregion(rgs, ind, edge)
SHEDLINE_ID = 0;
rgsId = rgs.Ids(ind);
rgsInd = rgs.Seg == rgsId;
rgsSpd = mean(rgs.SegValue(rgsInd));
rgsLocs = regionprops(rgsInd, 'pixellist');
% rgsXMin = min(rgsLocs.PixelList(:,2));
% rgsXMax = max(rgsLocs.PixelList(:,2));
% find upstream/downstream neighbours
% expCon = ones(3,1);
% expCon = [1 1 0; 0 1 0; 0 1 1];
expCon = ones(3); %[0 1 0; 1 1 1; 0 1 0];
rgsExpand = imdilate(rgsInd, expCon);
% delete extended pixels which are on edges
rgsExpand(edge) = 0;
% extend one more time
rgsExpand = imdilate(rgsExpand, expCon);
rgsExpand(edge) = 0;

rgsExpandIds = unique(rgs.Seg(rgsExpand));
rgsExpandIds(rgsExpandIds == rgsId) = []; % remove this region's id
rgsExpandIds(rgsExpandIds == SHEDLINE_ID) = []; % remove watershed line's id
nNeighbors = numel(rgsExpandIds);

% disp(nNeighbors)
if nNeighbors == 0
    return
end

% order these rgsExpandIds in increasing "closeness" - which is measured 
% based on the difference between a neighbour with the processed region.
spds = zeros(1, nNeighbors);
spdDif = zeros(1, nNeighbors);
for iId = 1:nNeighbors
    id = rgsExpandIds(iId);
    regionSpd = rgs.SegValue(rgs.Seg == id);
    spds(iId) = mean(regionSpd);
    spdDif(iId) = abs(spds(iId) - rgsSpd);
end
[~, spdOrder] = sort(spdDif, 'ascend');
rgsExpandIds = rgsExpandIds(spdOrder);

% % COMBINE ALL NEIGHBORS
% rgsMerge = ismember(rgs.Seg, [rgsExpandIds(:)', rgsId]);
% rgsMergeExpd = imdilate(rgsMerge, ones(3)); % ro remove boundaries
% rgsMerge = imerode(rgsMergeExpd, ones(3));
% imagesc(rgsMerge);

% allSeg = rgs.Seg ~= 120; allSeg = imdilate(allSeg, ones(3));
% otherRgs = ~ismember(rgs.Seg, [rgsExpandIds(:)', 120, rgsId]);
% otherRgs = imdilate(otherRgs, ones(3));
% rgsMerge = allSeg - otherRgs;
% rgsMerge(edge) = 0;

% COMBINE NEIGHBORS
isMerge = false;
rgsBase = imdilate(rgsInd, expCon); % dilate the base region for 1 pixel
for iNgbr = 1:nNeighbors
    rgsNgbr = rgs.Seg == rgsExpandIds(iNgbr);
    rgsNgbrExpd = imdilate(rgsNgbr, expCon);
    rgsMerge = imerode(rgsBase | rgsNgbrExpd, expCon);
    comBndry = rgsMerge & ~rgsInd & ~rgsNgbr;
    comBndryLen = sum(comBndry(:)); % number of pixel in the boundary
    comBndryOnEdge = sum(edge(comBndry));
    if comBndryOnEdge / comBndryLen < 0.3
        % merge
        % debug
%         subplot(2,2,1); imagesc(rgsInd);
%         subplot(2,2,2); imagesc(rgsNgbr);
%         subplot(2,2,3); imagesc(rgsMerge);
        %
        rgsMerge = rgsBase | rgsNgbrExpd;
        rgsBase = rgsMerge; % the core gets larger
        rgsInd = imerode(rgsMerge, expCon);
        rgs.Seg(logical(rgsInd)) = rgsId;
%         imagesc(rgsInd);
        rgs.isProcessed(rgs.Ids == rgsExpandIds(iNgbr)) = 1;
        isMerge = true;
    end
end

if ~isMerge
    return
end
% rgs.Seg(logical(rgsMerge)) = rgsId;
% rgs.isProcessed(ismember(rgs.Ids, rgsExpandIds)) = 1;
% continue merging
rgs = pr_growandmergeregion(rgs,ind, edge);
end
%%
function W = pr_cleansegmentation(W, M)
%pr_cleansegmentation decides whether to keep a segment
ids = unique(W(:));
ids(ids == 0) = [];
for id = ids'
    m = W == id;
    sin = sum(M(m));
    sall = sum(m(:));
    if 2*sin < sall
        W(m) = 0;
    end
end
end