import math

# Step 1: Define the properties of the two metric spaces.
# For the interval [0, 1], its length L is 1.
L = 1.0

# For the unit circle, its radius r is 1.
r = 1.0

# Step 2: Apply the known formula for the Gromov-Hausdorff distance
# between an interval of length L and a circle of radius r.
# The formula is d_GH = max(L/2, r).

# Calculate the term for the interval
half_L = L / 2

# Calculate the Gromov-Hausdorff distance
gh_distance = max(half_L, r)

# Step 3: Print the explanation and the final equation with the numbers.
print("The Gromov-Hausdorff distance between an interval of length L and a circle of radius r is given by the formula: max(L/2, r)")
print(f"For our case, L = {L} and r = {r}.")
print("The calculation is:")
print(f"max({half_L}, {r}) = {gh_distance}")
