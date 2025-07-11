import numpy as np

# --- Data from the problem ---
# n_carbons represents the number of carbon atoms in the R group of the ligand.
# distances represents the distance between adjacent polymer chains in Angstroms (Å).
n_carbons = np.array([1, 2, 3])
distances = np.array([12.4, 12.6, 13.2])
target_n_carbons = 8

# --- Step 1: State the problem and assumptions ---
print("The problem is to predict the inter-chain distance for a polymer formed with ligand 8 (n-octyl group).")
print("Ligand 8 has an alkyl chain with 8 carbon atoms.")
print("We assume that ligand 8 also forms a one-dimensional polymer, similar to ligands 1, 2, and 3.\n")

# --- Step 2: Analyze the trend in the provided data ---
print("--- Analyzing the Trend ---")
increase_1_to_2 = distances[1] - distances[0]
increase_2_to_3 = distances[2] - distances[1]
print(f"The distance increase from n=1 to n=2 carbons is: {distances[1]} - {distances[0]} = {increase_1_to_2:.1f} Å")
print(f"The distance increase from n=2 to n=3 carbons is: {distances[2]} - {distances[1]} = {increase_2_to_3:.1f} Å")
print("The rate of increase is accelerating.\n")

# --- Step 3: Create a model for extrapolation ---
print("--- Modeling for Extrapolation ---")
# Calculate the average slope using linear regression on all available data points.
avg_slope, _ = np.polyfit(n_carbons, distances, 1)
print(f"The average slope across the first three points (from linear regression) is {avg_slope:.2f} Å/carbon.")

# The most recent slope is the change from n=2 to n=3.
recent_slope = increase_2_to_3
print(f"The most recent slope (n=2 to n=3) is {recent_slope:.2f} Å/carbon.")

# A reasonable estimate for the ongoing slope is the average of the overall trend and the recent trend.
estimated_slope = (avg_slope + recent_slope) / 2
print(f"We'll use an estimated slope which is the average of these two values: ({avg_slope:.2f} + {recent_slope:.2f}) / 2 = {estimated_slope:.2f} Å/carbon.\n")

# --- Step 4: Predict the distance for n=8 ---
print("--- Predicting the Distance for Ligand 8 (n=8) ---")
# We extrapolate from the last known data point (n=3, distance=13.2 Å).
# The formula is: Predicted_Distance = Last_Distance + slope * (change in n)
last_n = n_carbons[-1]
last_distance = distances[-1]
change_in_n = target_n_carbons - last_n
predicted_distance = last_distance + estimated_slope * change_in_n

print(f"We predict the distance for n = {target_n_carbons} using the last data point (n={last_n}, d={last_distance} Å):")
print(f"Final Equation: D({target_n_carbons}) = D({last_n}) + slope * ({target_n_carbons} - {last_n})")
print(f"Calculation: {predicted_distance:.1f} = {last_distance} + {estimated_slope:.1f} * ({target_n_carbons} - {last_n})")
print(f"The predicted distance for the polymer with ligand 8 is {predicted_distance:.1f} Å.")
print("\nThis corresponds to a one-dimensional polymer with a distance of 15.7Å.")
