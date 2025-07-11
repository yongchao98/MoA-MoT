import numpy as np

# A list of the study sites for context
sites = [
    "Tena, Ecuador", "Cali, Colombia", "Barro Colorado Island, Panama",
    "Sarapiqui, Costa Rica", "Managua, Nicaragua", "Yoro, Honduras",
    "Belmopan, Belize", "Merida, Mexico", "Miami, Florida",
    "Charlotte, North Carolina", "State college, Pennsylvania", "Ottawa, Canada"
]

# Approximate latitudes (in degrees North) for the sites
latitudes = np.array([0, 3, 9, 10, 12, 15, 17, 21, 26, 35, 41, 45])

# --- Step 1: Define Initial State ---
# We'll create a hypothetical diversity dataset that models a strong latitudinal gradient.
# Let's assume a simple linear model: initial_diversity = 300 - 5 * latitude
initial_diversity = 300 - 5 * latitudes

print("--- Initial State Before Invasion ---")
print(f"Site Latitudes: {latitudes}")
print(f"Initial Alpha Diversity: {initial_diversity}")

# --- Step 2: Calculate Initial Slope ---
# We perform a linear regression (fitting a line of the form y = mx + c)
# to find the slope of diversity vs. latitude.
initial_slope, initial_intercept = np.polyfit(latitudes, initial_diversity, 1)

print("\nInitial linear equation describing diversity gradient:")
print(f"diversity = ({initial_slope:.2f} * latitude) + {initial_intercept:.2f}")

# --- Step 3: Model Post-Invasion State ---
# According to neutral theory, a successful invader causes biotic homogenization.
# The diversity at all sites drops and becomes similar.
# We model this by setting the final diversity to a constant low value, e.g., 20.
final_diversity = np.full_like(latitudes, 20)

print("\n--- Final State After Biotic Homogenization ---")
print(f"Site Latitudes: {latitudes}")
print(f"Final Alpha Diversity: {final_diversity}")

# --- Step 4: Calculate Final Slope ---
# We calculate the slope of the new, flattened diversity gradient.
final_slope, final_intercept = np.polyfit(latitudes, final_diversity, 1)

print("\nFinal linear equation describing diversity gradient:")
print(f"diversity = ({final_slope:.2f} * latitude) + {final_intercept:.2f}")

# --- Conclusion ---
print("\n--- Conclusion ---")
print("As predicted by Hubbell's unified theory, the successful invader causes biotic homogenization.")
print(f"The initial steep negative slope of {initial_slope:.2f} flattens, approaching {final_slope:.2f}.")
print("This indicates the disappearance of the latitudinal diversity gradient.")