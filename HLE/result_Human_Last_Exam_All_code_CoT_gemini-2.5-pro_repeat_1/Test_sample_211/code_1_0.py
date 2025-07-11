import numpy as np

# Plan:
# 1. Define the 12 research sites and their approximate latitudes.
# 2. Create a hypothetical "initial" alpha diversity dataset based on the principle of the 
#    latitudinal diversity gradient (higher diversity at lower latitudes).
# 3. Calculate and print the slope of the initial diversity gradient.
# 4. Simulate the biotic homogenization effect caused by the invasive species. This is done by
#    reducing diversity at all sites, with a larger absolute reduction at the more diverse
#    tropical sites. This flattens the gradient.
# 5. Calculate and print the new slope and the full equation after the invasion.
# 6. Conclude what happens to the slope based on the simulation.

# Site names and their approximate absolute latitudes
sites = {
    "Tena, Ecuador": 0.99,
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

latitudes = np.array(list(sites.values()))

# Step 2: Generate initial diversity data (higher diversity at lower latitudes)
# We'll use a simple linear model: Diversity = 250 - 4 * Latitude
initial_intercept = 250
initial_slope_param = -4
initial_diversity = initial_intercept + initial_slope_param * latitudes + np.random.normal(0, 5, len(latitudes)) # add some noise

# Step 3: Calculate the initial slope using linear regression
initial_slope, initial_b = np.polyfit(latitudes, initial_diversity, 1)

print("--- Initial State (Pre-Invasion) ---")
print(f"The initial slope of diversity vs. latitude is: {initial_slope:.4f}")
print("This negative slope indicates that insect diversity is higher closer to the equator (lower latitudes).")
print("\n")

# Step 4: Simulate the invasion's homogenizing effect
# The new diversity will be a fraction of the old diversity, plus a constant base.
# This reduces the overall spread and makes the gradient flatter.
# New Diversity = 0.5 * Old Diversity + 30
post_invasion_diversity = 0.5 * initial_diversity + 30

# Step 5: Calculate the new slope after the invasion
final_slope, final_intercept = np.polyfit(latitudes, post_invasion_diversity, 1)

print("--- Final State (Post-Invasion) ---")
print("After the invasive species causes biotic homogenization, the new relationship is calculated.")
print(f"The final slope of diversity vs. latitude is: {final_slope:.4f}")
print("The slope is now less steep (closer to zero), indicating a weaker latitudinal gradient.")
print("\nFinal Equation:")
# Output each number in the final equation
print(f"Diversity = {final_slope:.2f} * Latitude + {final_intercept:.2f}")

<<<The slope of insect diversity will decrease (become less steep and closer to zero).>>>