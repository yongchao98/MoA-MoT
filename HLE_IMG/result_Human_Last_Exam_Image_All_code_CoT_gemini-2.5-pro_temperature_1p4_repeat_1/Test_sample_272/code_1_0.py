import numpy as np

# Define the energy scales from the problem
# Current scale (Lambda_1) is 8 TeV
lambda_1_tev = 8.0
# Proposed new scale (Lambda_2) is 1.1 PeV
lambda_2_pev = 1.1

# Convert the scales to a common unit, GeV, for calculation
# 1 TeV = 1,000 GeV
# 1 PeV = 1,000,000 GeV
lambda_1_gev = lambda_1_tev * 1e3
lambda_2_gev = lambda_2_pev * 1e6

# As derived in the explanation, the fine-tuning measure Δ is directly proportional to the cutoff scale Λ.
# Therefore, the multiplicative factor by which Δ changes is the ratio of the new scale to the old scale.
# Factor = Δ_2 / Δ_1 = (k * Λ_2) / (k * Λ_1) = Λ_2 / Λ_1
multiplicative_factor = lambda_2_gev / lambda_1_gev

# --- Output the result ---
print("The problem asks for the multiplicative factor by which the fine-tuning measure changes.")
print("The fine-tuning measure Δ is proportional to the new physics cutoff scale Λ.")
print("The multiplicative factor is therefore the ratio of the new scale (Λ_2) to the old scale (Λ_1).\n")

print("The old scale is Λ_1 = 8 TeV = 8000 GeV.")
print("The new scale is Λ_2 = 1.1 PeV = 1100000 GeV.\n")

print("The calculation of the multiplicative factor is:")
print(f"Factor = Λ_2 / Λ_1 = {int(lambda_2_gev)} / {int(lambda_1_gev)}")

# Round the final answer to two decimal places as requested.
rounded_factor = round(multiplicative_factor, 2)

print(f"\nThe result is: {rounded_factor}")

print("\n<<<137.50>>>")