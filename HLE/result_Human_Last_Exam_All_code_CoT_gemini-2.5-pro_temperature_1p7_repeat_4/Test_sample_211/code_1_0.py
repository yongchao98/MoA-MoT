import numpy as np

# This script simulates the effect of an invasive species on the latitudinal diversity gradient
# according to Hubbell's unified neutral theory.

# The theory predicts biotic homogenization, where a successful invasive species
# reduces biodiversity more significantly in high-diversity regions (like the tropics)
# than in low-diversity regions. This causes the diversity gradient across latitudes
# to become less pronounced, or "flatten".

# Step 1: Define the research sites and their approximate latitudes.
# The absolute value of latitude is used as a measure of distance from the equator.
sites = {
    "Tena, Ecuador": -0.99,
    "Cali, Colombia": 3.45,
    "Barro Colorado Island, Panama": 9.15,
    "Sarapiqui, Costa Rica": 10.46,
    "Managua, Nicaragua": 12.13,
    "Yoro, Honduras": 15.14,
    "Belmopan, Belize": 17.25,
    "Merida, Mexico": 20.97,
    "Miami, Florida": 25.76,
    "Charlotte, North Carolina": 35.23,
    "State college, Pennsylvania": 40.80,
    "Ottawa, Canada": 45.42
}

latitudes = np.array([abs(lat) for lat in sites.values()])
site_names = list(sites.keys())

# Step 2: Model the initial diversity gradient.
# We'll assume a simple linear model where alpha diversity is highest at the equator
# and decreases with latitude.
# Model: Diversity = 500 - 10 * Latitude
initial_intercept = 500.0
initial_slope = -10.0
initial_diversity = initial_intercept + initial_slope * latitudes

# Step 3: Model the impact of the invasive species.
# The homogenization effect means diversity is reduced more in the tropics.
# We model the final diversity to show a greater reduction for higher initial values.
# Model: final_diversity = 0.5 * initial_diversity + 20
# Example (Tropics): Initial diversity of 500 becomes 0.5*500 + 20 = 270 (a drop of 230).
# Example (Temperate): Initial diversity of 50 becomes 0.5*50 + 20 = 45 (a drop of 5).
final_diversity = 0.5 * initial_diversity + 20

# Step 4: Calculate the slopes of the initial and final diversity gradients
# using linear regression to find the 'best-fit' line.
initial_fit_slope, initial_fit_intercept = np.polyfit(latitudes, initial_diversity, 1)
final_fit_slope, final_fit_intercept = np.polyfit(latitudes, final_diversity, 1)

# Step 5: Print the results and the conclusion.
print("Hubbell's unified neutral theory suggests that a widespread invasive species will cause biotic homogenization,")
print("ultimately weakening the latitudinal diversity gradient.\n")

print("--- Initial State ---")
print("The initial relationship between diversity and latitude is modeled as:")
print(f"Alpha Diversity = {initial_fit_intercept:.2f} + ({initial_fit_slope:.2f}) * Latitude")
print(f"The initial slope is {initial_fit_slope:.2f}\n")


print("--- Final State (after invasion) ---")
print("The final relationship after the invasion causes homogenization is:")
print(f"Alpha Diversity = {final_fit_intercept:.2f} + ({final_fit_slope:.2f}) * Latitude")
print(f"The final slope is {final_fit_slope:.2f}\n")

print("--- Conclusion ---")
print(f"The slope of diversity changed from {initial_fit_slope:.2f} to {final_fit_slope:.2f}.")
print("The magnitude of the slope has decreased, meaning the gradient has become flatter.")

<<<The slope of diversity will flatten, becoming less steep and approaching zero.>>>