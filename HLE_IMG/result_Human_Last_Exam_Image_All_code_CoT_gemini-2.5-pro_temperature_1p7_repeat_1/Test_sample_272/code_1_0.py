import math

# Define the old and new energy scales for new physics
lambda_old_tev = 8.0  # in TeV
lambda_new_pev = 1.1  # in PeV

# Convert the new scale to TeV to have consistent units
# 1 PeV = 1000 TeV
lambda_new_tev = lambda_new_pev * 1000

# The fine-tuning measure Delta is proportional to the cutoff scale Lambda.
# Therefore, the multiplicative factor of the change is the ratio of the new scale to the old scale.
# The other provided values (masses, coupling) are part of a proportionality constant that cancels out in the ratio.
factor = lambda_new_tev / lambda_old_tev

# Round the result to two decimal places as requested
factor_rounded = round(factor, 2)

# Print the final equation and the result
print(f"The old energy scale is {lambda_old_tev} TeV.")
print(f"The new proposed energy scale is {lambda_new_pev} PeV, which is {lambda_new_tev} TeV.")
print("The multiplicative factor is the ratio of these two scales.")
print(f"Final calculation: {lambda_new_tev} / {lambda_old_tev} = {factor_rounded}")
print(f"The fine-tuning measure changes by a multiplicative factor of {factor_rounded}.")

# Final answer in the required format
# <<<137.50>>>