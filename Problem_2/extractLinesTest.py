import ExtractLines
import numpy as np

#Two points, there must be a line between them

MIN_SEG_LENGTH = 0.05  # minimum length of each line segment (m)
LINE_POINT_DIST_THRESHOLD = 0.02  # max distance of pt from line to split
MIN_POINTS_PER_SEGMENT = 4  # minimum number of points per line segment
MAX_P2P_DIST = 1.0  # max distance between two adjent pts within a segment


params = {'MIN_SEG_LENGTH': MIN_SEG_LENGTH,
              'LINE_POINT_DIST_THRESHOLD': LINE_POINT_DIST_THRESHOLD,
              'MIN_POINTS_PER_SEGMENT': MIN_POINTS_PER_SEGMENT,
              'MAX_P2P_DIST': MAX_P2P_DIST}


test_theta = np.array([0., np.pi/8, np.pi/4])
test_rho = np.array([0.5, 2.3, 1.5])

alpha, r = ExtractLines.FitLine(test_theta, test_rho)

s = ExtractLines.FindSplit(test_theta, test_rho, alpha, r, params)

print ("Alpha: ", alpha)
print ("r: ", r)

print ("S: ", s)