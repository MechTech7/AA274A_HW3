import ExtractLines
import numpy as np

#Two points, there must be a line between them
test_theta = np.array([0., np.pi/4])
test_rho = np.array([0.5, 1.5])

alpha, r = ExtractLines.FitLine(test_theta, test_rho)

print ("Alpha: ", alpha)
print ("r: ", r)