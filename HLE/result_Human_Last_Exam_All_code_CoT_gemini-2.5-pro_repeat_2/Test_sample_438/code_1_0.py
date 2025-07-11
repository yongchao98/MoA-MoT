# Define the baseline relative risk and the fold increase.
# We can represent the baseline risk as 1 for simplicity.
baseline_relative_risk = 1

# The problem states the probability increases by 1000 fold.
fold_increase = 1000

# Calculate the new relative risk.
new_relative_risk = baseline_relative_risk * fold_increase

# Print the explanation and the final equation.
print("An HLA variant that increases the presentation of a disease-causing self-antigen would increase the risk of disease.")
print("This can be modeled as a multiplication of a baseline relative risk by the fold increase.")
print("\nCalculating the new relative risk:")
print(f"Baseline Relative Risk * Fold Increase = New Relative Risk")
print(f"{baseline_relative_risk} * {fold_increase} = {new_relative_risk}")