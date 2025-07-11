import numpy as np

# Define the initial and final energy scales.
# Initial scale Lambda_1 in TeV
lambda_1_tev = 8.0

# Final scale Lambda_2 in PeV
lambda_2_pev = 1.1

# Convert PeV to TeV for consistent units (1 PeV = 1000 TeV)
lambda_2_tev = lambda_2_pev * 1000

# As derived in the explanation, the multiplicative factor is the ratio of the two scales.
factor = lambda_2_tev / lambda_1_tev

# Print the calculation and the final result rounded to two decimal places.
print(f"The multiplicative factor is the ratio of the new scale to the old scale.")
print(f"Factor = {lambda_2_tev} TeV / {lambda_1_tev} TeV")
print(f"Result = {factor:.2f}")

# Final answer in the required format
print(f"\n<<< {factor:.2f} >>>")