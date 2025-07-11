# Baseline relative probability of antigen presentation. We can set this to 1 for a baseline reference.
baseline_presentation_probability = 1

# The fold increase caused by the HLA variant.
fold_increase = 1000

# Calculate the new relative probability of antigen presentation.
new_presentation_probability = baseline_presentation_probability * fold_increase

# Print the final equation showing the calculation.
print("The risk of immune activation is amplified. The calculation is:")
print(f"Final Equation: {baseline_presentation_probability} * {fold_increase} = {new_presentation_probability}")