import sys

# This script illustrates how a change in antigen presentation affects disease risk.

# Define the given factor of increase.
fold_increase = 1000

# Let's assume a baseline relative risk for an individual without this specific HLA variant.
# This is a hypothetical small number to represent a low starting risk.
base_relative_risk = 0.001 # or 0.1%

# For this model, we'll assume the change in disease risk is directly
# proportional to the change in the presentation of the disease-causing self-antigen.
# Calculate the new relative risk factor based on the fold increase.
new_relative_risk = base_relative_risk * fold_increase

# --- Output ---
print("This scenario describes how a genetic variation in an HLA gene affects disease risk.")
print("The presentation of a 'disease-causing self-antigen' is a key step in triggering autoimmunity.")
sys.stdout.flush()

print("\nLet's model the change in risk:")
print(f"A hypothetical baseline relative risk: {base_relative_risk}")
print(f"The fold increase in antigen presentation from the variant: {fold_increase}")

print("\nWe can calculate a new relative risk factor with the equation:")
# The final equation with each number is printed below as requested.
print(f"New Relative Risk = Baseline Relative Risk * Fold Increase")
print(f"{new_relative_risk:.1f} = {base_relative_risk} * {fold_increase}")

print("\nA 1000-fold increase in the presentation of a self-antigen that causes disease makes it far more likely")
print("that the immune system will recognize and attack the body's own tissues.")
print("\nConclusion: This would likely INCREASE a person's risk of developing the disease.")
sys.stdout.flush()
