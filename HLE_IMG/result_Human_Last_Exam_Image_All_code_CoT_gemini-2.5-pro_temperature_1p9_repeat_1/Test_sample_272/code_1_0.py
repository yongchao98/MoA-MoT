import numpy as np

# Step 1: Define the energy scales given in the problem.
lambda_1_tev = 8.0  # Current scale in TeV
lambda_2_pev = 1.1 # Proposed new scale in PeV

# Step 2: Convert the new scale from PeV to TeV to have consistent units.
# 1 PeV = 1000 TeV
lambda_2_tev = lambda_2_pev * 1000

# Step 3: As derived in the plan, the multiplicative factor is the ratio of the new scale to the old scale.
# This is because the fine-tuning measure Delta is directly proportional to the cutoff scale Lambda.
# Factor = Delta_2 / Delta_1 = Lambda_2 / Lambda_1
multiplicative_factor = lambda_2_tev / lambda_1_tev

# Step 4: Print the final calculation, showing each number in the equation,
# and round the result to two decimal places as requested.
print("The multiplicative factor is calculated as the ratio of the new scale to the old scale:")
print(f"Factor = {lambda_2_tev:.1f} TeV / {lambda_1_tev:.1f} TeV")
print(f"The result is: {multiplicative_factor:.2f}")
