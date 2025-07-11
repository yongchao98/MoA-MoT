import numpy as np

# Step 1: Define the known data points
# n = number of carbon atoms in the R group
n_values = np.array([1, 2, 3])
# d = distance between polymer chains in Angstroms
d_values = np.array([12.4, 12.6, 13.2])

print("Data points (n, d):")
for n, d in zip(n_values, d_values):
    print(f"n = {n}, d = {d} Å")

# Step 2: Perform a linear regression (fitting to a 1st-degree polynomial)
# This finds the slope (m) and intercept (c) for the line d = m*n + c
m, c = np.polyfit(n_values, d_values, 1)

print("\nAssuming a linear relationship d = m*n + c.")
print(f"The calculated slope (m) is: {m:.4f}")
print(f"The calculated intercept (c) is: {c:.4f}")
print(f"The best-fit linear equation is: d = {m:.4f}*n + {c:.4f}")

# Step 3: Predict the distance for ligand 8 (n=8)
n_predict = 8
d_predict = m * n_predict + c

print(f"\nPredicting the distance for ligand 8, where n = {n_predict}:")
# We need to print the equation with all the numbers.
# d(8) = m * 8 + c
print(f"d({n_predict}) = {m:.4f} * {n_predict} + {c:.4f}")
print(f"d({n_predict}) = {m * n_predict:.4f} + {c:.4f}")
print(f"d({n_predict}) = {d_predict:.4f} Å")


# Step 4: Compare with answer choices
# The question states that the product is a 1D polymer, and we must find the distance.
# Let's list the relevant answer choices:
# B. 14Å, C. 15.7Å, D. 17Å
print(f"\nThe predicted distance is approximately {d_predict:.1f} Å.")
print("Comparing this to the answer choices, the closest value is 15.7 Å (Choice C).")
