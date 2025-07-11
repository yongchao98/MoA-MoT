import numpy as np

# Define the 12 research sites and their approximate latitudes
sites = {
    "Tena, Ecuador": -0.99,
    "Cali, Colombia": 3.45,
    "Barro Colorado Island, Panama": 9.15,
    "Sarapiqui, Costa Rica": 10.47,
    "Managua, Nicaragua": 12.13,
    "Yoro, Honduras": 15.14,
    "Belmopan, Belize": 17.25,
    "Merida, Mexico": 20.96,
    "Miami, Florida": 25.76,
    "Charlotte, North Carolina": 35.22,
    "State College, Pennsylvania": 40.79,
    "Ottawa, Canada": 45.42
}

# Create arrays for latitude and a hypothetical initial alpha diversity
# Initial diversity shows a clear gradient, decreasing away from the equator
latitudes = np.array(list(sites.values()))
initial_diversity = np.array([150, 140, 125, 115, 110, 100, 95, 80, 60, 45, 30, 20])

# --- Step 1: Analyze the Initial State ---
# Calculate the slope of the initial diversity gradient using linear regression
initial_slope, initial_intercept = np.polyfit(latitudes, initial_diversity, 1)

print("--- Initial State ---")
print(f"Latitudes: {latitudes}")
print(f"Initial Diversities: {initial_diversity}")
print(f"Initial Slope of Diversity Gradient: {initial_slope:.4f}")
print("The initial negative slope shows that diversity is highest near the equator and decreases towards the poles.\n")


# --- Step 2: Model the Ultimate Outcome under Neutral Theory ---
# According to the theory, the single invasive species will ultimately out-compete
# all native species through stochastic drift, making the alpha diversity at every site equal to 1.
final_diversity = np.ones_like(latitudes)

# Calculate the slope of the final diversity "gradient"
final_slope, final_intercept = np.polyfit(latitudes, final_diversity, 1)

print("--- Final State (Ultimate Outcome) ---")
print("Under Hubbell's neutral theory, the invasive species drifts to fixation at all sites.")
print(f"Final Diversities: {final_diversity}")

# --- Step 3: Print the Final Equation and Conclusion ---
# The final equation is Diversity = (slope * Latitude) + intercept
print("\nThe final relationship between diversity and latitude can be described by the equation:")
print(f"Diversity = ({final_slope:.4f} * Latitude) + {final_intercept:.4f}")
print("Since the slope is effectively zero, the equation simplifies to: Diversity = 1")
print("\nConclusion: The slope of the diversity gradient will ultimately flatten and become zero.")

<<<0>>>