import numpy as np

# The problem provides data correlating the number of carbon atoms (n) in the ligand's alkyl chain
# with the distance (d) between adjacent polymer chains.
# We are given data for n=1, 2, and 3, and asked to predict the distance for n=8.

# Step 1: Define the known data points.
# n_values: Number of carbons in the alkyl chain.
# d_values: Corresponding distance in Angstroms (Å).
n_values = np.array([1, 2, 3])
d_values = np.array([12.4, 12.6, 13.2])

print("Analyzing the relationship between alkyl chain length and inter-polymer distance.")
print("Given data points (Number of Carbons, Distance in Å):")
for n, d in zip(n_values, d_values):
    print(f"({n}, {d})")

# Step 2: Model the data using a linear relationship, d = m*n + c.
# We use numpy's polyfit to perform a linear regression and find the slope (m) and intercept (c).
slope, intercept = np.polyfit(n_values, d_values, 1)

print("\nAssuming a linear model: distance = slope * (number of carbons) + intercept")
print(f"The calculated slope (m) is: {slope:.4f}")
print(f"The calculated intercept (c) is: {intercept:.4f}")

# Step 3: Predict the distance for ligand 8 (n-octyl), which has 8 carbon atoms.
n_predict = 8
d_predict = slope * n_predict + intercept

print(f"\nPredicting the distance for n = {n_predict}:")
# The final equation with all numbers plugged in
print(f"distance = {slope:.4f} * {n_predict} + {intercept:.4f}")
print(f"Predicted distance = {d_predict:.4f} Å")

# The product is a one-dimensional polymer with a distance between chains of approximately 15.13 Å.
# Comparing this value to the options, 15.7 Å (Choice C) is the closest.
print("\nThe calculated value of 15.133 Å is closest to the answer choice C (15.7 Å).")
