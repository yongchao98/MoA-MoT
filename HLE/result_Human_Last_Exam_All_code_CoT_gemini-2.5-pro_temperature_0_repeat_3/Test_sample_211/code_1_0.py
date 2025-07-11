import numpy as np

# 1. Define the sites and their approximate absolute latitudes
# Data from south to north
sites = [
    "Tena, Ecuador", "Cali, Colombia", "Barro Colorado Island, Panama",
    "Sarapiqui, Costa Rica", "Managua, Nicaragua", "Yoro, Honduras",
    "Belmopan, Belize", "Merida, Mexico", "Miami, Florida",
    "Charlotte, North Carolina", "State college, Pennsylvania", "Ottawa, Canada"
]
# Using absolute latitude as the gradient runs from the equator outwards
latitudes = np.array([
    0.99, 3.45, 9.15, 10.4, 12.13, 15.14, 17.25, 20.97, 25.76, 35.23, 40.80, 45.42
])

# 2. Assign a plausible initial alpha diversity to each site to create a gradient
# Higher diversity at lower latitudes
initial_diversity = np.array([
    250, 240, 220, 210, 200, 185, 170, 150, 120, 80, 65, 50
])

# 3. Calculate the initial slope of diversity vs. latitude
initial_slope, initial_intercept = np.polyfit(latitudes, initial_diversity, 1)

print("--- Initial State Before Invasion ---")
print(f"The initial slope of the diversity gradient is: {initial_slope:.2f}")
print("The equation describing the initial relationship is:")
print(f"Diversity = {initial_slope:.2f} * Latitude + {initial_intercept:.2f}")
print("-" * 35)

# 4. Simulate the invasion causing biotic homogenization
# We model this as a proportional loss of 60% of species at each site.
# Richer sites lose more species in absolute terms, flattening the gradient.
invasion_impact_factor = 0.60 # 60% of species are outcompeted
final_diversity = initial_diversity * (1 - invasion_impact_factor)

# 5. Calculate the final slope after the invasion
final_slope, final_intercept = np.polyfit(latitudes, final_diversity, 1)

print("\n--- Final State After Invasion ---")
print("The invasive species causes biotic homogenization, reducing diversity at all sites.")
print(f"The final slope of the diversity gradient is: {final_slope:.2f}")
print("The equation describing the final relationship is:")
# The prompt requires printing each number in the final equation
print(f"Diversity = {final_slope:.2f} * Latitude + {final_intercept:.2f}")
print("-" * 35)

# 6. Conclusion
print("\nConclusion:")
print(f"The slope changed from {initial_slope:.2f} to {final_slope:.2f}.")
print("As predicted by neutral theory, the slope has become less negative (flatter) and is closer to zero.")
print("This reflects the homogenization of the insect communities across the latitudinal gradient.")
