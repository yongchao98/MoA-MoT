import math

# Define the energy scales
# Current scale in TeV
lambda_1_tev = 8.0
# Proposed new scale in PeV
lambda_2_pev = 1.1

# Convert the new scale from PeV to TeV (1 PeV = 1000 TeV)
lambda_2_tev = lambda_2_pev * 1000

# The fine-tuning measure Delta is proportional to the bare Higgs mass, which is
# proportional to the cutoff scale Lambda. Therefore, the multiplicative factor
# is the ratio of the new scale to the old scale.
multiplicative_factor = lambda_2_tev / lambda_1_tev

# Print the explanation and the final calculation
print("The fine-tuning measure Δ is proportional to the cutoff scale Λ.")
print("The multiplicative change is the ratio of the new scale to the old scale.")
print(f"Factor = Λ_new / Λ_old")
print(f"Factor = {lambda_2_tev} TeV / {lambda_1_tev} TeV")
print(f"The multiplicative factor is: {multiplicative_factor:.2f}")
