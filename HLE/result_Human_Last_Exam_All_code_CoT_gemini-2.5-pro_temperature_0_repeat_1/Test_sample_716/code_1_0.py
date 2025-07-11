import math

# The curve is defined by the parametric equations:
# x(t) = sin(pi * t)
# y(t) = sin(t)
# z(t) = cos(2 * t)

# The Hausdorff dimension of a smooth curve is equal to its topological dimension, which is 1.
# We verify the curve is smooth by checking its derivative.
# The derivative vector is r'(t) = (pi*cos(pi*t), cos(t), -2*sin(2*t)).
# Since all components of the derivative are continuous functions, the curve is smooth.
# Therefore, its Hausdorff dimension is 1.

hausdorff_dimension = 1

# The final equation is simply that the dimension is 1.
# We print the number from this conclusion.
print(f"The final equation is: Dimension = {hausdorff_dimension}")
print(f"The Hausdorff dimension of the curve is: {hausdorff_dimension}")
