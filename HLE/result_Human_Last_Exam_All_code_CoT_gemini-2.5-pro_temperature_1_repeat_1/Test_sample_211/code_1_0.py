import numpy as np

# Step 1: Define the sites and their approximate latitudes.
sites = {
    "Tena, Ecuador": -0.99,
    "Cali, Colombia": 3.45,
    "Barro Colorado Island, Panama": 9.15,
    "Sarapiqui, Costa Rica": 10.4,
    "Managua, Nicaragua": 12.13,
    "Yoro, Honduras": 15.14,
    "Belmopan, Belize": 17.25,
    "Merida, Mexico": 20.97,
    "Miami, Florida": 25.76,
    "Charlotte, North Carolina": 35.23,
    "State College, Pennsylvania": 40.80,
    "Ottawa, Canada": 45.42
}

# Use the absolute value of latitude for the gradient calculation
latitudes = np.array([abs(lat) for lat in sites.values()])
site_names = list(sites.keys())

# Step 2: Create a hypothetical "pre-invasion" alpha diversity dataset.
# We model a strong latitudinal gradient: Diversity = 250 - 4 * Latitude
# This creates a steep negative slope.
pre_invasion_intercept = 250
pre_invasion_slope_factor = -4
pre_invasion_diversity = pre_invasion_intercept + pre_invasion_slope_factor * latitudes

# Step 3: Model the "post-invasion" state based on neutral theory.
# Biotic homogenization flattens the slope. We'll model this by
# reducing the slope's magnitude significantly. The overall diversity also drops.
# New Diversity = 180 - 1.5 * Latitude
post_invasion_intercept = 180
post_invasion_slope_factor = -1.5
post_invasion_diversity = post_invasion_intercept + post_invasion_slope_factor * latitudes

# Step 4: Calculate the slope and intercept from the data using linear regression.
# This confirms the slopes of our models.
pre_slope, pre_intercept = np.polyfit(latitudes, pre_invasion_diversity, 1)
post_slope, post_intercept = np.polyfit(latitudes, post_invasion_diversity, 1)

# Step 5: Print the results, including the equation for each scenario.
print("--- Pre-Invasion Scenario ---")
print("This scenario represents a strong Latitudinal Diversity Gradient.")
print("The calculated linear equation for diversity is:")
print(f"Alpha Diversity = ({pre_slope:.2f} * Latitude) + {pre_intercept:.2f}")
print("\n--- Post-Invasion Scenario (Ultimate Outcome) ---")
print("This scenario represents the effect of biotic homogenization under neutral theory.")
print("The calculated linear equation for diversity is:")
print(f"Alpha Diversity = ({post_slope:.2f} * Latitude) + {post_intercept:.2f}")
print("\n--- Conclusion ---")
print(f"The slope of diversity changed from {pre_slope:.2f} to {post_slope:.2f}.")
print("The magnitude of the slope has decreased, meaning the gradient has become flatter and less steep.")
