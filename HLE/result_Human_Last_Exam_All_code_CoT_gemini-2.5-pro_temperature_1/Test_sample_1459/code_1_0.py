import math

# The Gromov-Hausdorff distance between the interval [0, 1] and the unit circle
# with the intrinsic metric is given by the formula (pi - 1) / 2.

# diam_interval is the diameter of the interval [0, 1]
diam_interval = 1

# diam_circle is the diameter of the unit circle (radius 1) with the intrinsic metric
diam_circle = math.pi

# The Gromov-Hausdorff distance is calculated based on the analysis.
# The lower bound is 0.5 * |diam_circle - diam_interval|.
# An upper bound is found to be the same value by constructing a specific map.
distance = (diam_circle - diam_interval) / 2

# We print the final equation with the numbers plugged in.
print(f"The Gromov-Hausdorff distance is given by the equation:")
print(f"d_GH = (π - 1) / 2")
print(f"d_GH = ({diam_circle} - {diam_interval}) / 2")
print(f"d_GH = {distance}")
