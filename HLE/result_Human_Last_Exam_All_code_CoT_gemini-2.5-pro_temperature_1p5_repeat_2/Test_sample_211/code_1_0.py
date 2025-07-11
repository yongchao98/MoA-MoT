import numpy as np

# Step 1: Define the study sites and their approximate absolute latitudes.
# The gradient is a function of distance from the equator (latitude 0).
sites_and_latitudes = {
    "Tena, Ecuador": 0.99,
    "Cali, Colombia": 3.45,
    "Barro Colorado Island, Panama": 9.15,
    "Sarapiqui, Costa Rica": 10.46,
    "Managua, Nicaragua": 12.13,
    "Yoro, Honduras": 15.14,
    "Belmopan, Belize": 17.25,
    "Merida, Mexico": 20.96,
    "Miami, Florida": 25.76,
    "Charlotte, North Carolina": 35.22,
    "State college, Pennsylvania": 40.79,
    "Ottawa, Canada": 45.42
}
latitudes = np.array(list(sites_and_latitudes.values()))

# Step 2: Simulate the initial state with a strong latitudinal diversity gradient.
# We model diversity as a linear function of latitude: Diversity = Intercept + Slope * Latitude
# A large negative slope represents a steep gradient. We add some noise for realism.
np.random.seed(42)  # for reproducible results
initial_intercept = 180
initial_slope = -3.0
initial_noise = np.random.normal(0, 8, len(latitudes))
initial_diversity = initial_intercept + initial_slope * latitudes + initial_noise

# Calculate the best-fit line for the initial data
initial_slope_fit, initial_intercept_fit = np.polyfit(latitudes, initial_diversity, 1)

print("--- Initial State Before Invasive Species ---")
print("A strong latitudinal diversity gradient is observed.")
print(f"The equation for the initial diversity gradient is: Diversity = {initial_intercept_fit:.2f} + ({initial_slope_fit:.2f}) * Latitude")
print("\n")


# Step 3: Simulate the final state after the invasive species causes biotic homogenization.
# The gradient flattens, meaning the slope becomes less negative (closer to zero).
# Overall diversity is also reduced.
final_intercept = 90
final_slope = -0.5
final_noise = np.random.normal(0, 8, len(latitudes))
final_diversity = final_intercept + final_slope * latitudes + final_noise

# Calculate the best-fit line for the final data
final_slope_fit, final_intercept_fit = np.polyfit(latitudes, final_diversity, 1)

print("--- Final State After Invasive Species Spreads ---")
print("Biotic homogenization has occurred, making communities more similar.")
print(f"The equation for the final diversity gradient is: Diversity = {final_intercept_fit:.2f} + ({final_slope_fit:.2f}) * Latitude")
print("\n")

# Step 4: Output the conclusion by comparing the slopes.
print("--- Conclusion ---")
print(f"The slope of the latitudinal diversity gradient has changed from {initial_slope_fit:.2f} to {final_slope_fit:.2f}.")
print("The magnitude of the slope has decreased significantly, indicating that the diversity gradient has flattened.")
