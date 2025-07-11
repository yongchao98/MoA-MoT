import numpy as np
from scipy import stats

# 1. Define the study sites and their approximate latitudes (in degrees).
sites = {
    "Tena, Ecuador": 0.99,
    "Cali, Colombia": 3.44,
    "Barro Colorado Island, Panama": 9.15,
    "Sarapiqui, Costa Rica": 10.42,
    "Managua, Nicaragua": 12.13,
    "Yoro, Honduras": 15.14,
    "Belmopan, Belize": 17.25,
    "Merida, Mexico": 20.97,
    "Miami, Florida": 25.76,
    "Charlotte, North Carolina": 35.22,
    "State College, Pennsylvania": 40.79,
    "Ottawa, Canada": 45.42
}
latitudes = np.array(list(sites.values()))
site_names = list(sites.keys())

# 2. Generate initial alpha diversity values to model the latitudinal gradient.
# We'll create a strong negative correlation with latitude, where diversity is highest
# at the equator (low latitude) and decreases northward.
# y = intercept + slope * x
initial_intercept_theoretical = 250
initial_slope_theoretical = -4.5
# Add some random noise to make the data more realistic.
np.random.seed(42)
noise = np.random.normal(0, 10, len(latitudes))
initial_diversity = initial_intercept_theoretical + initial_slope_theoretical * latitudes + noise
initial_diversity[initial_diversity < 15] = 15  # Ensure diversity is positive.

# 3. Simulate the effect of the invasive species.
# Under neutral theory, the invader causes biotic homogenization, reducing diversity
# across all sites and making them more similar. We model this by shrinking the
# diversity values towards a new, lower mean.
final_diversity = initial_diversity * 0.3 + 25

# 4. Calculate the slope of the diversity gradient before and after the invasion.
initial_slope, initial_intercept, _, _, _ = stats.linregress(latitudes, initial_diversity)
final_slope, final_intercept, _, _, _ = stats.linregress(latitudes, final_diversity)

# 5. Print the explanation and results.
print("Prediction based on Hubbell's Unified Neutral Theory")
print("=" * 60)
print("Hubbell's theory posits that species are ecologically equivalent. A successful invader spreads and, through random demographic processes, displaces native species.")
print("This leads to 'biotic homogenization,' reducing overall diversity and decreasing the differences between communities.")
print("\nThis simulation demonstrates the effect on the latitudinal diversity gradient:")

print("\n--- Initial State ---")
print("Initially, there is a strong latitudinal gradient in diversity.")
print(f"The equation for the initial diversity gradient is:")
print(f"Insect Diversity = {initial_intercept:.2f} + ({initial_slope:.2f}) * Latitude")


print("\n--- Post-Invasion State ---")
print("After the invasion, diversity is reduced and homogenized across all sites.")
print(f"The equation for the final diversity gradient is:")
print(f"Insect Diversity = {final_intercept:.2f} + ({final_slope:.2f}) * Latitude")

print("\n--- Conclusion ---")
print(f"The slope of diversity versus latitude has changed from {initial_slope:.2f} to {final_slope:.2f}.")
print("The magnitude of the slope has decreased significantly, bringing it closer to zero.")
print("Therefore, under Hubbell's theory, the invasive species will cause the slope of insect diversity to flatten.")
print("=" * 60)
<<<The slope will flatten, approaching zero.>>>