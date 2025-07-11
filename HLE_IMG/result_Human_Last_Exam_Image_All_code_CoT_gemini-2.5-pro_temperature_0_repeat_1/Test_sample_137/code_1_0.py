import numpy as np
from scipy import stats

# Step 1: Define the known data from the problem description.
# The R groups for ligands 1, 2, and 3 are Methyl, Ethyl, and n-Propyl.
# The number of carbon atoms in these groups are 1, 2, and 3, respectively.
num_carbons = np.array([1, 2, 3])

# The corresponding distances between adjacent polymer chains are given in Angstroms (Å).
distances = np.array([12.4, 12.6, 13.2])

print("--- Analysis of Polymer Chain Distance vs. Alkyl Chain Length ---")
print("Known data points (Number of Carbons, Distance in Å):")
print(f"Ligand 1 (1 Carbon): {distances[0]} Å")
print(f"Ligand 2 (2 Carbons): {distances[1]} Å")
print(f"Ligand 3 (3 Carbons): {distances[2]} Å")

# Step 2: Assume a linear relationship and perform a linear regression to find the trend.
# The model is: distance = slope * num_carbons + intercept
slope, intercept, r_value, p_value, std_err = stats.linregress(num_carbons, distances)

print("\nStep 2: Modeling the trend with a best-fit line.")
print(f"The derived linear model is: distance = {slope:.3f} * (num_carbons) + {intercept:.3f}")

# Step 3: Use the model to predict the distance for ligand 8.
# Ligand 8 has an n-octyl group (R = n-octyl), which has 8 carbon atoms.
carbons_to_predict = 8
predicted_distance = slope * carbons_to_predict + intercept

print(f"\nStep 3: Predicting the distance for Ligand 8 ({carbons_to_predict} carbons).")
# As requested, showing the final equation with all numbers.
print("Final Equation:")
print(f"Predicted Distance = {slope:.3f} * {carbons_to_predict} + {intercept:.3f}")
print(f"Result: {predicted_distance:.3f} Å")

print("\n--- Conclusion ---")
print(f"The predicted distance is approximately {predicted_distance:.2f} Å.")
print("The problem states that analogous one-dimensional polymers are formed.")
print("Comparing our prediction to the answer choices, the closest value is 15.7 Å (Choice C).")
print("Therefore, the expected product is a one-dimensional polymer with a distance of 15.7 Å between chains.")
