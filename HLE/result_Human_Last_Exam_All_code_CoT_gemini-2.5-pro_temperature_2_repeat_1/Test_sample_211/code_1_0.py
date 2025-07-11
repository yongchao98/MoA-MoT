import numpy as np
from scipy import stats

# 1. Model the Initial State
# List of sites and their approximate latitudes (in degrees north)
sites = {
    "Tena, Ecuador": -0.5,
    "Cali, Colombia": 3.4,
    "Barro Colorado Island, Panama": 9.1,
    "Sarapiqui, Costa Rica": 10.4,
    "Managua, Nicaragua": 12.1,
    "Yoro, Honduras": 15.1,
    "Belmopan, Belize": 17.2,
    "Merida, Mexico": 21.0,
    "Miami, Florida": 25.8,
    "Charlotte, North Carolina": 35.2,
    "State College, Pennsylvania": 40.8,
    "Ottawa, Canada": 45.4
}

latitudes = np.array(list(sites.values()))

# Generate initial diversity with a strong latitudinal gradient
# Diversity = Intercept - Slope * Latitude. Let's use Intercept=200, Slope=3.5
initial_slope_param = -3.5
initial_intercept_param = 200
# Add a small amount of random noise to make it more realistic
np.random.seed(0)
noise = np.random.normal(0, 5, len(latitudes))
initial_diversity = initial_intercept_param + initial_slope_param * latitudes + noise

# 2. Calculate the Initial Slope
initial_slope, initial_intercept, _, _, _ = stats.linregress(latitudes, initial_diversity)

print("--- Initial State (Pre-Invasion) ---")
print("The relationship between alpha diversity and latitude is modeled by the equation:")
print(f"Diversity = {initial_slope:.2f} * Latitude + {initial_intercept:.2f}")
print("\nThe initial slope shows a strong negative relationship, with diversity decreasing as latitude increases.\n")

# 3. Model the Post-Invasion State
# The invasive species causes biotic homogenization, driving diversity down at all sites
# to a similar, low level, largely independent of latitude.
# Let's model the new diversity as being randomly distributed around a low value (e.g., 15).
np.random.seed(42)
final_diversity = np.random.normal(15, 3, len(latitudes))

# 4. Calculate the Final Slope
final_slope, final_intercept, _, _, _ = stats.linregress(latitudes, final_diversity)

print("--- Final State (Post-Invasion) ---")
print("After homogenization by the invasive species, the new relationship is:")
print(f"Diversity = {final_slope:.2f} * Latitude + {final_intercept:.2f}")
print("\nThe final slope is much closer to zero, indicating the latitudinal gradient has been erased.")
