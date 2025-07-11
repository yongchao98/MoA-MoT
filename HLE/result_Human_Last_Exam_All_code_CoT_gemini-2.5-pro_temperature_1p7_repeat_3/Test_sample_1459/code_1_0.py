import math

# 1. Define the properties of the two metric spaces.
# Space X: The interval [0, 1] with the standard absolute value metric.
# The diameter of X is the maximum distance between its endpoints.
diam_x = 1.0

# Space Y: The unit circle with the intrinsic metric (shortest arc length).
# A unit circle has radius 1. Its circumference is 2 * pi.
# The diameter of Y is the distance between antipodal points, which is
# half the circumference.
diam_y = math.pi

# 2. Use the formula for the lower bound of the Gromov-Hausdorff distance.
# d_GH(X, Y) >= |diam(X) - diam(Y)| / 2
# For the given spaces, this bound is tight and gives the exact distance.
distance = abs(diam_x - diam_y) / 2

# 3. Print the step-by-step calculation and the final equation.
print("This script calculates the Gromov-Hausdorff distance d_GH([0, 1], S^1).")
print("-" * 70)
print(f"Step 1: Determine the diameter of each space.")
print(f"  - For the interval [0, 1], the diameter is {diam_x}.")
print(f"  - For the unit circle (radius 1), the circumference is 2*pi.")
print(f"    The diameter is half the circumference: pi ≈ {diam_y:.5f}.")

print("\nStep 2: Apply the formula for the Gromov-Hausdorff distance.")
print("  d_GH(X, Y) = |Diameter(X) - Diameter(Y)| / 2")

print("\nStep 3: Substitute the values and compute.")
print(f"  d_GH([0, 1], S^1) = |{diam_x} - {diam_y:.5f}| / 2")
print(f"                  = {abs(diam_x - diam_y):.5f} / 2")
print(f"                  = {distance:.5f}")

print("\nStep 4: The final equation with symbolic pi.")
# In the final equation, we show each number involved.
number_1 = 1
number_2 = 2
print(f"d_GH([0, 1], S^1) = (pi - {number_1}) / {number_2}")