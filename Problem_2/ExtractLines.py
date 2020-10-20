#!/usr/bin/env python

############################################################
# ExtractLines.py
#
# This script reads in range data from a csv file, and
# implements a split-and-merge to extract meaningful lines
# in the environment.
############################################################

# Imports
import numpy as np
from PlotFunctions import *


############################################################
# functions
############################################################

def ExtractLines(RangeData, params):
    '''
    This function implements a split-and-merge line extraction algorithm.

    Inputs:
        RangeData: (x_r, y_r, theta, rho)
            x_r: robot's x position (m).
            y_r: robot's y position (m).
            theta: (1D) np array of angle 'theta' from data (rads).
            rho: (1D) np array of distance 'rho' from data (m).
        params: dictionary of parameters for line extraction.
    Outputs:
        alpha: (1D) np array of 'alpha' for each fitted line (rads).
        r: (1D) np array of 'r' for each fitted line (m).
        segend: np array (N_lines, 4) of line segment endpoints. Each row represents [x1, y1, x2, y2].
        pointIdx: (N_lines,2) segment's first and last point index.
    '''

    # Extract useful variables from RangeData
    x_r = RangeData[0]
    y_r = RangeData[1]
    theta = RangeData[2]
    rho = RangeData[3]

    ### Split Lines ###
    N_pts = len(rho)
    r = np.zeros(0)
    alpha = np.zeros(0)
    pointIdx = np.zeros((0, 2), dtype=np.int)

    # This implementation pre-prepartitions the data according to the "MAX_P2P_DIST"
    # parameter. It forces line segmentation at sufficiently large range jumps.
    rho_diff = np.abs(rho[1:] - rho[:(len(rho)-1)])
    LineBreak = np.hstack((np.where(rho_diff > params['MAX_P2P_DIST'])[0]+1, N_pts))
    startIdx = 0
    for endIdx in LineBreak:
        alpha_seg, r_seg, pointIdx_seg = SplitLinesRecursive(theta, rho, startIdx, endIdx, params)
        N_lines = r_seg.size

        ### Merge Lines ###
        if (N_lines > 1):
            alpha_seg, r_seg, pointIdx_seg = MergeColinearNeigbors(theta, rho, alpha_seg, r_seg, pointIdx_seg, params)
        r = np.append(r, r_seg)
        alpha = np.append(alpha, alpha_seg)
        pointIdx = np.vstack((pointIdx, pointIdx_seg))
        startIdx = endIdx

    N_lines = alpha.size

    ### Compute endpoints/lengths of the segments ###
    segend = np.zeros((N_lines, 4))
    seglen = np.zeros(N_lines)
    for i in range(N_lines):
        rho1 = r[i]/np.cos(theta[pointIdx[i, 0]]-alpha[i])
        rho2 = r[i]/np.cos(theta[pointIdx[i, 1]-1]-alpha[i])
        x1 = rho1*np.cos(theta[pointIdx[i, 0]])
        y1 = rho1*np.sin(theta[pointIdx[i, 0]])
        x2 = rho2*np.cos(theta[pointIdx[i, 1]-1])
        y2 = rho2*np.sin(theta[pointIdx[i, 1]-1])
        segend[i, :] = np.hstack((x1, y1, x2, y2))
        seglen[i] = np.linalg.norm(segend[i, 0:2] - segend[i, 2:4])

    ### Filter Lines ###
    # Find and remove line segments that are too short
    goodSegIdx = np.where((seglen >= params['MIN_SEG_LENGTH']) &
                          (pointIdx[:, 1] - pointIdx[:, 0] >= params['MIN_POINTS_PER_SEGMENT']))[0]
    pointIdx = pointIdx[goodSegIdx, :]
    alpha = alpha[goodSegIdx]
    r = r[goodSegIdx]
    segend = segend[goodSegIdx, :]

    # change back to scene coordinates
    segend[:, (0, 2)] = segend[:, (0, 2)] + x_r
    segend[:, (1, 3)] = segend[:, (1, 3)] + y_r

    return alpha, r, segend, pointIdx

#params["MIN_POINTS_PER_SEGMENT"]

'''def SplitLinesRecursive(theta, rho, startIdx, endIdx, params):
    
    This function executes a recursive line-slitting algorithm, which
    recursively sub-divides line segments until no further splitting is
    required.

    Inputs:
        theta: (1D) np array of angle 'theta' from data (rads).
        rho: (1D) np array of distance 'rho' from data (m).
        startIdx: starting index of segment to be split.
        endIdx: ending index of segment to be split.
        params: dictionary of parameters.
    Outputs:
        alpha: (1D) np array of 'alpha' for each fitted line (rads).
        r: (1D) np array of 'r' for each fitted line (m).
        idx: (N_lines,2) segment's first and last point index.

    HINT: Call FitLine() to fit individual line segments.
    HINT: Call FindSplit() to find an index to split at.
    
    ########## Code starts here ##########
    alpha, r = FitLine(theta[startIdx:endIdx+1], rho[startIdx:endIdx+1])

    if (endIdx - startIdx) + 1 < params["MIN_POINTS_PER_SEGMENT"]:
        return alpha, r, np.array([[startIdx, endIdx]])

    splitIdx = FindSplit(theta[startIdx:endIdx+1], rho[startIdx:endIdx+1], alpha, r, params)
    if splitIdx < 0:
        return alpha, r, np.array([[startIdx, endIdx]])
    
    alpha_1, r_1, idx_1 = SplitLinesRecursive(theta, rho, startIdx, startIdx + splitIdx, params)
    alpha_2, r_2, idx_2 = SplitLinesRecursive(theta, rho, startIdx + splitIdx, endIdx, params)

    idx = np.append(idx_1, idx_2, axis=0)
    alpha = np.append(alpha_1, alpha_2)
    r = np.append(r_1, r_2)
    ########## Code ends here ##########

    
    return alpha, r, idx'''
'''def SplitLinesRecursive(theta, rho, startIdx, endIdx, params):
    alpha, r = FitLine(theta[startIdx:endIdx], rho[startIdx:endIdx])

    if (endIdx - startIdx) < params["MIN_POINTS_PER_SEGMENT"]:
        return np.array([alpha]), np.array([r]), np.array([[startIdx, endIdx]])
    
    s = FindSplit(theta[startIdx:endIdx], rho[startIdx:endIdx], alpha, r, params)

    if s == -1:
        return np.array([alpha]), np.array([r]), np.array([[startIdx, endIdx]])
    
    alpha_1, r_1, idx_1 = SplitLinesRecursive(theta, rho, startIdx, startIdx + s, params)
    alpha_2, r_2, idx_2 = SplitLinesRecursive(theta, rho, startIdx + s, endIdx, params)

    alpha_op = np.hstack((alpha_1, alpha_2))
    r_op = np.hstack((r_1, r_2))
    idx_op = np.vstack((idx_1, idx_2))

    print (alpha_op.shape)
    return alpha_op, r_op, idx_op'''

def SplitLinesRecursive(theta, rho, startIdx, endIdx, params):
    theta_line = theta[startIdx:endIdx]
    rho_line = rho[startIdx:endIdx]

    alpha, r = FitLine(theta_line, rho_line)
    
    if (endIdx-startIdx) <= params["MIN_POINTS_PER_SEGMENT"]:
        return np.array([alpha]), np.array([r]), np.array([[startIdx, endIdx]])

    s = FindSplit(theta_line, rho_line, alpha, r, params)

    if s == -1:
        print ("Cutoff because of me")
        return np.array([alpha]), np.array([r]), np.array([[startIdx, endIdx]])

    alpha_1, r_1, idx_1 = SplitLinesRecursive(theta, rho, startIdx, startIdx + s, params)
    alpha_2, r_2, idx_2 = SplitLinesRecursive(theta, rho, startIdx + s, endIdx, params)

    alpha_out = np.hstack((alpha_1, alpha_2))
    r_out = np.hstack((r_1, r_2))
    idx_out = np.vstack((idx_1, idx_2))

    return alpha_out, r_out, idx_out

def FindSplit(theta, rho, alpha, r, params):
    '''
    This function takes in a line segment and outputs the best index at which to
    split the segment, or -1 if no split should be made.

    The best point to split at is the one whose distance from the line is
    the farthest, as long as this maximum distance exceeds
    LINE_POINT_DIST_THRESHOLD and also does not divide the line into segments
    smaller than MIN_POINTS_PER_SEGMENT. Otherwise, no split should be made.

    Inputs:
        theta: (1D) np array of angle 'theta' from data (rads).
        rho: (1D) np array of distance 'rho' from data (m).
        alpha: 'alpha' of input line segment (1 number).
        r: 'r' of input line segment (1 number).
        params: dictionary of parameters.
    Output:
        splitIdx: idx at which to split line (return -1 if it cannot be split).
    '''
    ########## Code starts here ##########
    cos_subt = np.cos(theta - alpha)
    dist = np.multiply(rho, cos_subt) - r 
    dist = np.abs(dist)

    length = dist.shape[0]
    min_points = params["MIN_POINTS_PER_SEGMENT"]
    dist[0:min_points] = 0.0
    dist[length - (min_points-1):] = 0.0
    
    splitIdx = np.argmax(dist)

    if dist[splitIdx] < params["LINE_POINT_DIST_THRESHOLD"]:
        return -1
    
    ########## Code ends here ##########
    return splitIdx

def FitLine(theta, rho):
    
    '''This function outputs a least squares best fit line to a segment of range
    data, expressed in polar form (alpha, r).

    Inputs:
        theta: (1D) np array of angle 'theta' from data (rads).
        rho: (1D) np array of distance 'rho' from data (m).
    Outputs:
        alpha: 'alpha' of best fit for range data (1 number) (rads).
        r: 'r' of best fit for range data (1 number) (m).'''
    
    ########## Code starts here ##########
    length = rho.shape[0]
    ##alpha Calculation:
    f_sum = np.sum(np.multiply(np.square(rho), np.sin(2*theta)), axis=0)
    sec_sum = 0
    for i in range(length):
        for j in range(length):
            sec_sum += rho[i] * rho[j] * np.cos(theta[i])*np.sin(theta[j])
    
    sec_sum *= -2.0/length
    numerator = f_sum + sec_sum

    d_sum = np.sum(np.multiply(np.square(rho), np.cos(2*theta)), axis=0)
    d_sec_sum = 0
    for i in range(length):
        for j in range(length):
            d_sec_sum += rho[i] * rho[j] * np.cos(theta[i] + theta[j])
    
    d_sec_sum *= -1.0/length
    denominator = d_sum + d_sec_sum

    alpha = 0.5 * np.arctan2(numerator, denominator) + np.pi/2
    r_array = np.multiply(rho, np.cos(theta - alpha))

    r = 1.0/length * np.sum(r_array)
    ########## Code ends here ##########
    return alpha, r



def MergeColinearNeigbors(theta, rho, alpha, r, pointIdx, params):
    #TODO: Debug this, its very faulty
    '''
    This function merges neighboring segments that are colinear and outputs a
    new set of line segments.

    Inputs:
        theta: (1D) np array of angle 'theta' from data (rads).
        rho: (1D) np array of distance 'rho' from data (m).
        alpha: (1D) np array of 'alpha' for each fitted line (rads).
        r: (1D) np array of 'r' for each fitted line (m).
        pointIdx: (N_lines,2) segment's first and last point indices.
        params: dictionary of parameters.
    Outputs:
        alphaOut: output 'alpha' of merged lines (rads).
        rOut: output 'r' of merged lines (m).
        pointIdxOut: output start and end indices of merged line segments.

    HINT: loop through line segments and try to fit a line to data points from
          two adjacent segments. If this line cannot be split, then accept the
          merge. If it can be split, do not merge.
    '''
    ########## Code starts here ##########
    alphaOut = []
    rOut = []
    pointIdxOut = []

    length = alpha.shape[0]

    i = 0
    while i < length - 1:
        curr_idx = pointIdx[i]
        next_idx = pointIdx[i+1]

        theta_line = theta[curr_idx[0]:next_idx[1]]
        rho_line = rho[curr_idx[0]:next_idx[1]]

        t_alpha, t_r = FitLine(theta_line, rho_line)
        s = FindSplit(theta_line, rho_line, t_alpha, t_r, params)

        if s < 0:
            alphaOut.append(t_alpha)
            rOut.append(t_r)
            pointIdxOut.append(np.array([curr_idx[0], next_idx[1]]))
            i += 2
        else:
            alphaOut.append(alpha[i])
            rOut.append(r[i])
            pointIdxOut.append(pointIdx[i])
            i+=1

    alphaOut = np.array(alphaOut)
    rOut = np.array(rOut)
    pointIdxOut = np.array(pointIdxOut)


    ########## Code ends here ##########
    return alphaOut, rOut, pointIdxOut


#----------------------------------
# ImportRangeData
def ImportRangeData(filename):

    data = np.genfromtxt('./RangeData/'+filename, delimiter=',')
    x_r = data[0, 0]
    y_r = data[0, 1]
    theta = data[1:, 0]
    rho = data[1:, 1]
    return (x_r, y_r, theta, rho)
#----------------------------------


############################################################
# Main
############################################################
def main():
    # parameters for line extraction (mess with these!)
    MIN_SEG_LENGTH = 0.05  # minimum length of each line segment (m)
    LINE_POINT_DIST_THRESHOLD = 0.02  # max distance of pt from line to split
    MIN_POINTS_PER_SEGMENT = 4  # minimum number of points per line segment
    MAX_P2P_DIST = 0.5  # max distance between two adjent pts within a segment

    # Data files are formated as 'rangeData_<x_r>_<y_r>_N_pts.csv
    # where x_r is the robot's x position
    #       y_r is the robot's y position
    #       N_pts is the number of beams (e.g. 180 -> beams are 2deg apart)

    filename = 'rangeData_5_5_180.csv'
    # filename = 'rangeData_4_9_360.csv'
    # filename = 'rangeData_7_2_90.csv'

    # Import Range Data
    RangeData = ImportRangeData(filename)

    params = {'MIN_SEG_LENGTH': MIN_SEG_LENGTH,
              'LINE_POINT_DIST_THRESHOLD': LINE_POINT_DIST_THRESHOLD,
              'MIN_POINTS_PER_SEGMENT': MIN_POINTS_PER_SEGMENT,
              'MAX_P2P_DIST': MAX_P2P_DIST}

    #theta = RangeData[2]
    #rho = RangeData[3]

    #alpha, r = FitLine(theta, rho)

    alpha, r, segend, pointIdx = ExtractLines(RangeData, params)

    ax = PlotScene()
    ax = PlotData(RangeData, ax)
    ax = PlotRays(RangeData, ax)
    ax = PlotLines(segend, ax)

    plt.show(ax)

############################################################

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        import traceback
        traceback.print_exc()
