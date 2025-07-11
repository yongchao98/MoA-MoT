import numpy as np

# The curve is given by r(t) = (sin(pi*t), sin(t), cos(2*t)).
# To find the Hausdorff dimension, we first check if the curve is smooth.
# We do this by computing its derivative (the tangent vector) and checking if it's ever the zero vector.
# r'(t) = (pi*cos(pi*t), cos(t), -2*sin(2*t))
#
# For the tangent vector to be zero, all components must be zero:
# 1) pi * cos(pi * t) = 0  => t = k + 1/2 for integer k.
# 2) cos(t) = 0             => t = m*pi + pi/2 for integer m.
#
# For both to be true simultaneously, we would need k + 1/2 = m*pi + pi/2 for integers k and m.
# This requires k = (m + 1/2)*pi - 1/2.
# Since pi is irrational, this equation has no integer solutions for k and m.
#
# This means the tangent vector is never the zero vector. The curve is a regular, smooth curve.
# A smooth curve is a 1-dimensional manifold.
# The Hausdorff dimension of a k-dimensional smooth manifold is k.
# In our case, k=1.
hausdorff_dimension = 1

# The "final equation" is simply that the dimension is 1.
# The following code prints the numerical result of this analysis.
print(hausdorff_dimension)