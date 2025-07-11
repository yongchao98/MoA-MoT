import numpy as np

# Step 1: Define the study sites and their approximate latitudes (in degrees).
sites = {
    "Tena, Ecuador": 1,
    "Cali, Colombia": 3,
    "Barro Colorado Island, Panama": 9,
    "Sarapiqui, Costa Rica": 10,
    "Managua, Nicaragua": 12,
    "Yoro, Honduras": 15,
    "Belmopan, Belize": 17,
    "Merida, Mexico": 21,
    "Miami, Florida": 26,
    "Charlotte, North Carolina": 35,
    "State College, Pennsylvania": 41,
    "Ottawa, Canada": 45
}

latitudes = np.array(list(sites.values()))

# Step 2: Create a hypothetical initial state with a strong latitudinal diversity gradient.
# Let's assume diversity decreases significantly with latitude.
# Initial Diversity = 250 - 4 * latitude (with some noise)
np.random.seed(0) # for reproducibility
initial_diversity = 250 - 4 * latitudes + np.random.normal(0, 5, len(latitudes))

# Step 3: Model the final state after the invasion leads to biotic homogenization.
# Overall diversity decreases, and the gradient becomes flatter (less steep).
# Final Diversity = 200 - 2 * latitude (with some noise)
final_diversity = 200 - 2 * latitudes + np.random.normal(0, 5, len(latitudes))

# Step 4: Calculate the slope and intercept for both states using linear regression.
initial_slope, initial_intercept = np.polyfit(latitudes, initial_diversity, 1)
final_slope, final_intercept = np.polyfit(latitudes, final_diversity, 1)

# Step 5: Print the results, showing the full equation for both states.
print("Initial State before homogenization:")
print(f"The equation for the diversity gradient is: Diversity = {initial_intercept:.2f} + ({initial_slope:.2f} * Latitude)")
print(f"The initial slope is {initial_slope:.2f}\n")

print("Final State after homogenization:")
print(f"The equation for the diversity gradient is: Diversity = {final_intercept:.2f} + ({final_slope:.2f} * Latitude)")
print(f"The final slope is {final_slope:.2f}\n")

print("Conclusion:")
print(f"The slope changes from {initial_slope:.2f} to {final_slope:.2f}.")
print("As predicted by neutral theory, the biotic homogenization caused by the invader makes the latitudinal diversity gradient less steep.")