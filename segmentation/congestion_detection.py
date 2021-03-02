"""
 @author Tin T. Nguyen
 @email t.t.nguyen-3@tudelft.nl
 @create date 2020-11-03
 @modify date 2020-11-05
 @desc [description]
""" 

import numpy as np
import pandas as pd
import warnings
from skimage.segmentation import chan_vese
from segmentation.chanvese.chanvese import chanvese
from skimage.morphology import label as skimage_label
from utils.preprocess import preprocess_speed

def detect_congestion(speed, method='Chan-Vese', version='shawnlankton', speed_threshold=50, debug=False):
    congestion_mask = None
    if method == 'Chan-Vese':
        un = speed / np.max(speed)
        inimask = np.int8(speed < speed_threshold)
        if not np.any(inimask):
            return np.zeros(speed.shape)
        if version == 'shawnlankton':
            congestion_mask, _, _ = chanvese(speed, init_mask=speed<speed_threshold, alpha=0, max_its=200)
        else:
            congestion_mask = chan_vese(un, mu=0, lambda1=1, lambda2=1, tol=1e-3, max_iter=500,dt=0.5, init_level_set=inimask, extended_output=False)
    return congestion_mask

def mark_congested_link(speed, 
                        link_route_df,
                        x,
                        x_raw,
                        congestion_detection_method='Chan-Vese',
                        version='shawnlankton',
                        sparse_time=60,
                        time_step=2):

    congestion_mask = np.zeros(speed.shape)
    # replace NaN (if any) with the maximum speed
    speed = preprocess_speed(speed)
    if speed is None:
        return congestion_mask
    # spatial resolution (m)
    dx = x[2] - x[1]
    # downsample
    speed = speed[:,0::time_step]
    if congestion_detection_method == 'naive':
        pass
    else:
        congestion_mask = detect_congestion(speed)
    if not np.any(congestion_mask):
        raise Exception('no congestion')
    
    def process_link(link):
        # upstream/downstream detectors
        detector_ind = np.where(np.logical_and((x_raw >= link['distance_along_route'][0]),
                                               (x_raw < link['distance_along_route'][-1])))[0]
        if np.any(detector_ind):
            # upstream
            det_ind_up = detector_ind[0]
            # the corresponding asm index
            asm_ind_up = min(np.round(x_raw[det_ind_up] / dx).astype(np.int), speed.shape[0]-1)
            # downstream
            det_ind_down = detector_ind[-1]
            # the corresponding asm index
            asm_ind_down = min(np.round(x_raw[det_ind_down] / dx).astype(np.int), speed.shape[0]-1)
        else:
            # no detectors on this link
            # use our best bet: asm data
            asm_ind_up = min(np.round(link['distance_along_route'][0]/dx).astype(np.int), speed.shape[0]-1)
            asm_ind_down = min(np.round(link['distance_along_route'][-1]/dx).astype(np.int), speed.shape[0]-1)
        # link congestion mask using detector locations
        link_congestion_matrix = congestion_mask[asm_ind_up:asm_ind_down+1,:]
        labeled_congestion_region = skimage_label(link_congestion_matrix)
        # upstream front
        link_congestion_upstream = labeled_congestion_region[0,:]
        # downstream front
        link_congestion_downstream = labeled_congestion_region[-1,:]
        # check if congestion involves the whole link
        # first, exaggerate upstream/downstream front to the whole link
        link_asm_ind = np.where(np.logical_and((x >= link['distance_along_route'][0]),
                                               (x < link['distance_along_route'][-1])))[0]
        all_relevant_ind = np.append(link_asm_ind, [asm_ind_up, asm_ind_down])
        link_asm_ind_up = np.min(all_relevant_ind)
        link_asm_ind_down = np.max(all_relevant_ind)
        link_congestion_mask = congestion_mask[link_asm_ind_up:link_asm_ind_down+1,:]
        link_congestion_any = np.any(link_congestion_mask, axis=0)
#         link_congestion_all = np.all(link_congestion_mask, axis=0)
        link_congestion_all = np.zeros(link_congestion_any.shape)
        L, nL = skimage_label(link_congestion_mask, return_num=True)
        for i in range(nL):
            cm = L == (i+1)
            if np.all(np.any(cm,axis=1)):
                link_congestion_all = np.logical_or(link_congestion_all, 
                                                    np.any(cm,axis=0))
        return pd.Series({'link_congestion_all': link_congestion_all, 
                         'link_congestion_any': link_congestion_any,
                         'link_congestion_upstream': link_congestion_upstream, 
                         'link_congestion_downstream': link_congestion_downstream})
    link_congestion_mark = link_route_df.apply(lambda l: process_link(l), axis=1)
    return link_congestion_mark, congestion_mask

# function m = cpc_hp_estimatecongestedmask(speed, method, debug)
# %CPC_HP_ESTIMATECONGESTEDMASK function binarizes traffic states into free-flow or congested.
# % Available methods: 
# %       bimo        : assumes traffic speeds follow a bimodal distribution
# %         lbimo     : location-wise
# %         gbimo     : roadstretch-wise
# %       adap        : estimates free speeds at each location
# %       contours    : active contours without edges (Chen-Vese)
# %%
# % m: congested mask
# switch(method)
#     case 'adap'
#         m = false(size(speed));
#         for iLoc = 1:size(speed,1)
#             [~, ~, loDev] = pr_estimatefree(speed(iLoc, :));
#             m(iLoc, :) = speed(iLoc,:) < loDev;
#         end
#         m = imopen(m, ones(11,1));
#         m = imopen(m, ones(1, 11));
#     case 'lbimo'
#         % bimodal at each location
#         m = false(size(speed));
#         for iLoc = 1:size(speed,1)
#             u = speed(iLoc,:);
#             nu = u ./ max(u);
#             m(iLoc, :) = ~imbinarize(nu);
#         end
#     case 'gbimo'
#         nu = speed ./ (max(speed(:)));
#         m = ~imbinarize(nu);
#     case 'gadapbimo'
#         % bimodal at each road stretch, which has the same free speed at
#         % all locations
#         % divide the road into multiple sections according to free speeds
#         [~, mask, nsection, range] = cp_estimatefreespeeds(speed);
#         m = zeros(size(speed));
#         for i = 1:nsection
#             % section locs
#             sx = range(i,1):range(i,2);
#             % section speeds
#             sspeed = speed(sx, :);
#             % normalised section speeds
#             nsu = sspeed ./ max(sspeed(:));
#             % binarise
#             sm = ~imbinarize(nsu);
#             % add to the overall pattern mask
#             m(sx,:) = sm;
#         end
#     case 'contour_wu'
#         un = speed / max(speed(:));
#         inimask = speed < 50;
#         m = chenvese(un,inimask,500,0,'chan', debug);
#     case 'contour_builtin'
#         un = speed / max(speed(:));
#         inimask = speed < 50;
#         m = activecontour(un,inimask,500,'Chan-Vese');
#         % m = pr_refinecongmask(speed, m);
#     otherwise
#         nu = speed ./ (max(speed(:)));
#         m = ~imbinarize(nu);
# end
# end

# function [u0, o, loDev] = pr_estimatefree(u)
# % u: speeds
# % p: prominences
# [pks,~,~,p] = findpeaks(u);
# pks = ceil(pks);
# p = ceil(p);

# uMax = max(u(:)); uMax = ceil(uMax);
# nP = length(pks);
# m = zeros(uMax, nP);
# for i = 1:nP
#     pV = pks(i); % peak's value
#     pO = max(p(i),1);   % peak's prominence
#     m(pV:-1:pV-pO+1, i) = 1;
# end
# o = sum(m,2); % overlapping
# % find the last peak
# o(end+1) = 0; % just in case the last value is a peak
# %% version 1
# % [~, l, ~, ~] = findpeaks(o, 'MinPeakHeight', 5);
# % u0 = l(end);
# % % deviation
# % % find all (speed) peaks which together with their offsets include u0
# % isRelPkInds = m(u0, :) ~= 0;
# % relM = m(:, isRelPkInds);
# % r = sum(relM, 2);
# % % loDev = find(r,1); % lower deviation
# % loDev = 2*u0 - max(pks);
# %% version 2

# [pc, l, ~, ~] = findpeaks(o);
# % find the peak at such sum of all the subsequent peaks is at least 5
# pc_r = flipud(pc); % reverse 'counting peaks' in pc
# l_r = flipud(l);
# c = cumsum(pc_r);
# i = find(c>=5, 1); % find the first one that is at least 5
# u0 = l_r(i);
# % deviation
# % find all (speed) peaks which together with their offsets include u0
# % isRelPkInds = m(u0, :) ~= 0;
# % relM = m(:, isRelPkInds);
# % r = sum(relM, 2);
# % loDev = find(r,1); % lower deviation
# loDev = 2*u0 - max(pks);
# end

# %% refine congestion mask by incorporating Watershed
# function refMask = pr_refinecongmask(speed, iniMask)
# % smooth speed to remove noise
# spdBlur = imgaussfilt(speed, 0.5);
# % gradient
# [Gmag, ~] = imgradient(spdBlur);
# % normalize Gmag
# Gmag = Gmag / max(Gmag(:));
# [spdEdge] = edge(speed, 'canny', [], 0.7);
# Gmag(spdEdge == 1) = 1;
# % Gmag = imhmin(Gmag, 0.015);
# EXTREME_CONGESTED_SPEED = 20;
# spdExtrCong = speed < EXTREME_CONGESTED_SPEED;
# Gmag(spdExtrCong) = min(Gmag(spdExtrCong));
# % spdFree = spdROI >= 100; %FREE_SPEED;
# % spdFree = imerode(spdFree, ones(5));
# Gmag(spdEdge == 1) = 1;
# % Gmag(spdFree) = min(Gmag(:));
# shrunkenFreeMask = imerode(~iniMask, ones(15));
# Gmag(shrunkenFreeMask) = min(Gmag(:));
# % Watershed
# % WGmag = watershed(Gmag, 4);
# WGmag = watershed(Gmag);
# %% remove free traffic regions
# % freeMask = data_hp_getpatternmask(speed, freeMask);
# % freeMask = ~spdFree; % this is temporary
# % WGmagFreeIds = unique(WGmag(~freeMask));
# % WGmagFreeIds = unique(WGmag(shrunkenFreeMask));
# % WGmag(ismember(WGmag, WGmagFreeIds)) = 0;
# % WGmag = double(WGmag) .* freeMask;
# freeSegIds = WGmag(shrunkenFreeMask);
# WGmag(ismember(WGmag, freeSegIds)) = 0;
# %% retain only congested regions
# congSegIds = WGmag(iniMask);
# refMask = ismember(WGmag, congSegIds);
# % remove shed lines
# refMask = imclose(refMask, ones(3));
# end