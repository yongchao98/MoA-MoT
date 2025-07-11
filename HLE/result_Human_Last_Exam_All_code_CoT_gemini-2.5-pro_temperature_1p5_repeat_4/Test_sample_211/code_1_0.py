import numpy as np

# Step 1: Define the study sites and their approximate latitudes.
# Southern latitudes are negative. We will use the absolute value for the gradient calculation
# as the gradient runs from the equator outwards.
sites = {
    "Tena, Ecuador": -0.99,
    "Cali, Colombia": 3.44,
    "Barro Colorado Island, Panama": 9.15,
    "Sarapiqui, Costa Rica": 10.42,
    "Managua, Nicaragua": 12.13,
    "Yoro, Honduras": 15.1,
    "Belmopan, Belize": 17.25,
    "Merida, Mexico": 20.97,
    "Miami, Florida": 25.76,
    "Charlotte, North Carolina": 35.23,
    "State college, Pennsylvania": 40.79,
    "Ottawa, Canada": 45.42
}

latitudes = np.array(list(sites.values()))

# Step 2: Simulate the initial state - a strong latitudinal diversity gradient.
# We create a hypothetical linear relationship where diversity decreases with latitude.
# Initial Diversity = 300 - 5 * |Latitude|
initial_diversity = 300 - 5 * np.abs(latitudes)

# Step 3: Calculate and display the initial slope of diversity vs. latitude.
# np.polyfit(x, y, 1) returns [slope, intercept] for a linear fit.
initial_slope, initial_intercept = np.polyfit(np.abs(latitudes), initial_diversity, 1)

print("--- Initial State (Pre-Invasion) ---")
print("A strong latitudinal diversity gradient is observed.")
print("The relationship between alpha diversity and absolute latitude is modeled by the equation:")
print(f"Diversity = {initial_slope:.2f} * |Latitude| + {initial_intercept:.2f}\n")


# Step 4: Simulate the effect of the invasive species causing biotic homogenization.
# Homogenization reduces diversity everywhere and makes communities more similar.
# This shrinks the diversity values towards a common, lower mean.
# This transformation will reduce high-diversity sites more than low-diversity sites, flattening the slope.
final_diversity = initial_diversity * 0.15 + 30

# Step 5: Calculate and display the final slope after homogenization.
final_slope, final_intercept = np.polyfit(np.abs(latitudes), final_diversity, 1)

print("--- Final State (Post-Invasion) ---")
print("After the invasive species causes biotic homogenization:")
print("The relationship between alpha diversity and absolute latitude has flattened. The new equation is:")
print(f"Diversity = {final_slope:.2f} * |Latitude| + {final_intercept:.2f}\n")


# Step 6: Explicitly state the conclusion.
print("--- Conclusion ---")
print(f"The initial slope's magnitude was {abs(initial_slope):.2f}.")
print(f"The final slope's magnitude is {abs(final_slope):.2f}.")
print("The slope has become much flatter, meaning its value is closer to zero.")
print("This demonstrates that the differences in diversity across the latitudinal gradient have been severely reduced.")
