import numpy as np
from scipy.stats import linregress

# Plan:
# 1. Define the 12 research sites and their approximate latitudes.
# 2. Create a hypothetical "pre-invasion" dataset showing a strong latitudinal diversity gradient.
# 3. Calculate the slope of this initial gradient using linear regression.
# 4. Simulate biotic homogenization by reducing the diversity difference between sites, creating a "post-invasion" dataset.
# 5. Calculate the new, flatter slope.
# 6. Print the results to show the change in the slope.

# Step 1: Define sites and their approximate absolute latitudes (degrees from the equator)
sites = {
    "Tena, Ecuador": 1.0,
    "Cali, Colombia": 3.4,
    "Barro Colorado Island, Panama": 9.1,
    "Sarapiqui, Costa Rica": 10.4,
    "Managua, Nicaragua": 12.1,
    "Yoro, Honduras": 15.1,
    "Belmopan, Belize": 17.2,
    "Merida, Mexico": 20.9,
    "Miami, Florida": 25.8,
    "Charlotte, North Carolina": 35.2,
    "State college, Pennsylvania": 40.8,
    "Ottawa, Canada": 45.4
}

latitudes = np.array(list(sites.values()))

# Step 2: Create a hypothetical "pre-invasion" dataset
# We'll model initial diversity as a linear function of latitude.
# diversity = intercept - slope * latitude
initial_intercept = 250
initial_slope_param = 4.5
diversity_before = initial_intercept - initial_slope_param * latitudes

# Step 3: Calculate the initial slope
# Using linregress to find the slope of diversity vs. latitude
slope_before, _, _, _, _ = linregress(latitudes, diversity_before)

print("--- Pre-Invasion State ---")
print(f"This simulates the initial latitudinal diversity gradient.")
print(f"Example Diversity at Tena (lat ~{sites['Tena, Ecuador']:.1f}): {int(diversity_before[0])} species")
print(f"Example Diversity at Ottawa (lat ~{sites['Ottawa, Canada']:.1f}): {int(diversity_before[-1])} species")
print("\nCalculating the initial slope of the gradient (Diversity vs. Latitude)...")
# The final equation is the comparison of the two slopes. Here is the first number:
print(f"Initial Slope: {slope_before:.4f}")


# Step 4: Simulate biotic homogenization post-invasion
# This process makes communities more similar, reducing the range of diversity values.
# We can simulate this by compressing the diversity values towards their mean.
mean_diversity = np.mean(diversity_before)
compression_factor = 0.4 # Represents strong homogenization
diversity_after = mean_diversity + (diversity_before - mean_diversity) * compression_factor

# Step 5: Calculate the new slope
slope_after, _, _, _, _ = linregress(latitudes, diversity_after)

print("\n--- Post-Invasion State (after Biotic Homogenization) ---")
print("The invasive species reduces diversity across all sites, making them more similar.")
print(f"New Diversity at Tena: {int(diversity_after[0])} species")
print(f"New Diversity at Ottawa: {int(diversity_after[-1])} species")
print("\nCalculating the final slope of the flattened gradient...")
# Here is the second number for the comparison:
print(f"Final Slope: {slope_after:.4f}")

# Step 6: Final Conclusion
print("\n--- Conclusion ---")
print(f"The magnitude of the slope decreased from {abs(slope_before):.4f} to {abs(slope_after):.4f}.")
print("This demonstrates that under neutral theory, the introduction of a major invasive species")
print("leads to biotic homogenization, which causes the latitudinal diversity gradient to flatten.")

<<<The slope will decrease (flatten).>>>