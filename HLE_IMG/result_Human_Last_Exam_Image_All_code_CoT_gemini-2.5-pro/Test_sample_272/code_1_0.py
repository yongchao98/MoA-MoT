import math

# Define the energy scales from the problem description
# The initial scale where no new physics has been found.
lambda_1_tev = 8.0

# The proposed new scale for a future machine.
lambda_2_pev = 1.1

# To calculate the ratio, we need the scales in the same units.
# We convert PeV to TeV using the relation: 1 PeV = 1000 TeV.
lambda_2_tev = lambda_2_pev * 1000

# The problem asks for the multiplicative factor by which the fine-tuning measure changes.
# As explained in the plan, the fine-tuning measure Delta is directly proportional
# to the cutoff scale Lambda. Therefore, the multiplicative factor is the ratio
# of the new scale to the old scale.
multiplicative_factor = lambda_2_tev / lambda_1_tev

# Print the final result, showing the equation as requested.
print(f"The calculation for the multiplicative factor is based on the ratio of the energy scales.")
print(f"Old scale: {lambda_1_tev} TeV")
print(f"New scale: {lambda_2_tev} TeV")
print(f"Multiplicative Factor = {lambda_2_tev} / {lambda_1_tev}")
print(f"The factor is: {multiplicative_factor:.2f}")