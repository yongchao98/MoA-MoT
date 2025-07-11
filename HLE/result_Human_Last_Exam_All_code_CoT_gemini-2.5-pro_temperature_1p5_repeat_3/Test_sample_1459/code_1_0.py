import math

# Step 1: Define the properties of the two spaces.
# Space X: Interval [0, 1]
# The metric is d(a,b) = |a-b|.
diam_X = 1.0

# Space Y: Unit circle (radius=1) with the intrinsic metric (arc length).
# The circumference is 2 * pi * r.
# The diameter is the distance between opposite points, which is half the circumference.
radius_Y = 1.0
circumference_Y = 2 * math.pi * radius_Y
diam_Y = circumference_Y / 2

# Step 2: Apply the formula for the Gromov-Hausdorff distance.
# For an interval and a circle, the distance is given by the diameter inequality,
# which in this case is not just a bound but the exact value.
# d_GH(X, Y) = |diam(X) - diam(Y)| / 2
gh_distance = abs(diam_X - diam_Y) / 2

# Step 3: Print the results step-by-step.
print(f"Let X be the interval [0, 1] and Y be the unit circle.")
print(f"The diameter of X is diam(X) = {diam_X}")
print(f"The diameter of Y is diam(Y) = pi * r = {diam_Y:.4f}")
print("\nThe Gromov-Hausdorff distance is calculated using the formula: |diam(X) - diam(Y)| / 2")
print(f"d_GH(X, Y) = |{diam_X} - {diam_Y:.4f}| / 2")
print(f"d_GH(X, Y) = {abs(diam_X - diam_Y):.4f} / 2")
print(f"d_GH(X, Y) = {gh_distance:.4f}")

# Final numerical answer to be captured
# print(f"<<<{gh_distance}>>>")