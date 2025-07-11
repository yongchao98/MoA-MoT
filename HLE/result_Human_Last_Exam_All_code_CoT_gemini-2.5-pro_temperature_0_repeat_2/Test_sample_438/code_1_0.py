# Define the baseline probability or risk factor. We can set it to 1 to represent a standard baseline.
initial_risk_factor = 1

# Define the fold increase caused by the HLA variant.
fold_increase = 1000

# Calculate the new, relative risk factor.
new_risk_factor = initial_risk_factor * fold_increase

# Print the final equation showing each number.
print(f"The calculation for the new relative risk is:")
print(f"{initial_risk_factor} (Baseline Risk) * {fold_increase} (Fold Increase) = {new_risk_factor} (New Relative Risk)")