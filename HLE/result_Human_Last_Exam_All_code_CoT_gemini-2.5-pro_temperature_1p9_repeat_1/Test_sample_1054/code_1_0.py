# NOTE: The data used in this script is a hypothetical example based on typical medical
# study values, as the actual data from your source was not provided.
# You can replace the data in the 'simulated_data' list with your own values.

# Step 1: Hypothesize a dataset.
# Each entry is a tuple: (Aortomesenteric Diameter cutoff in mm, Sensitivity in %, Specificity in %)
simulated_data = [
    (10.5, 55, 96),
    (9.8, 62, 92),
    (9.1, 75, 89),
    (8.4, 88, 81),
    (7.7, 94, 78),
    (7.0, 98, 65)
]

# Step 2: Define the required thresholds from the user's request.
sensitivity_threshold = 60
specificity_threshold = 80

print("Evaluating each AM diameter against the criteria:")
print(f"Required: Sensitivity > {sensitivity_threshold}% AND Specificity > {specificity_threshold}%\n")

# Step 3 & 4: Iterate through the data, check conditions, and collect valid diameters.
valid_diameters = []
for diameter, sensitivity, specificity in simulated_data:
    # Check if the data point meets both conditions
    if sensitivity > sensitivity_threshold and specificity > specificity_threshold:
        valid_diameters.append(diameter)
        print(f"Checking Diameter {diameter} mm -> (Sensitivity={sensitivity}%, Specificity={specificity}%) -> PASSED")
    else:
        print(f"Checking Diameter {diameter} mm -> (Sensitivity={sensitivity}%, Specificity={specificity}%) -> FAILED")

# Step 5 & 6: Determine the maximum value and print the final answer.
if valid_diameters:
    # The question "At most, what..." asks for the maximum qualifying value.
    max_valid_diameter = max(valid_diameters)
    print("\n-------------------------------------------------")
    print(f"The diameters that satisfy the conditions are: {valid_diameters}")
    print(f"The maximum (at most) Aortomesenteric diameter that meets the criteria is: {max_valid_diameter} mm")
    print("-------------------------------------------------")
else:
    print("\n-------------------------------------------------")
    print("No Aortomesenteric diameter in the dataset satisfies both conditions.")
    print("-------------------------------------------------")