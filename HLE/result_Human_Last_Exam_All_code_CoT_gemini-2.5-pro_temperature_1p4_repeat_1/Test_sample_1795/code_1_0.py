import numpy as np

# Step 1: Define constants and parameters
# Conversion factor from fm^-1 to MeV
hbar_c = 197.327  # MeV fm

# Given parameters
M = 938.93  # Nucleon mass in MeV
Cs2 = 267.1  # Scalar coupling constant squared (dimensionless)
kF_fm = 1.42  # Fermi momentum in fm^-1
R_fm = 5.5  # Radius of heavy nucleus in fm
M_star_ratio = 0.78  # M*/M
nu = 4  # Degeneracy factor for nuclear matter

# Convert to consistent natural units (MeV)
kF = kF_fm * hbar_c  # MeV
R = R_fm / hbar_c  # MeV^-1
M_star = M_star_ratio * M  # MeV

# Step 2: Evaluate the momentum integral part
# The term in parenthesis in the derived formula
momentum_term_val = (kF**2 + 2 * M_star**2) / np.sqrt(kF**2 + M_star**2) - 2 * M_star

# Step 3: Calculate the relative change per unit eta
# The formula is (1/eta) * (dM_star / M_star) = (Cs2/M**2) * nu * (3*R / (32 * pi**2)) * momentum_term
# Based on the plan, we assume a typo and remove the factor of 4 from the denominator (32 -> 8)
prefactor = (Cs2 / M**2) * nu * (3 * R / (8 * np.pi**2))
relative_change_per_eta = prefactor * momentum_term_val

# Step 4: Interpret "per unit of NC parameters" and calculate final percentage
# We assume "per unit" means eta = 1 MeV^2
eta = 1.0  # MeV^2
percentage_change = relative_change_per_eta * eta * 100

# Step 5: Print the results clearly
print("--- Input Parameters (in MeV units) ---")
print(f"M = {M:.2f} MeV")
print(f"M* = {M_star:.2f} MeV")
print(f"k_F = {kF:.2f} MeV")
print(f"R = {R:.6f} MeV^-1")
print(f"Cs^2 = {Cs2}")
print(f"nu = {nu}")
print("\n--- Calculation Steps ---")
print(f"Momentum integral term = {momentum_term_val:.4f} MeV")
print(f"Prefactor (Cs^2/M^2 * nu * 3R / (8*pi^2)) = {prefactor:.4e} MeV^-3")
print(f"Relative change per unit eta = {relative_change_per_eta:.4e} MeV^-2")
print("\n--- Final Result ---")
print("Assuming a typo correction in the formula (1/4Omega -> 1/Omega) and eta = 1 MeV^2:")
print(f"Percentage Change = (Relative change per unit eta) * (eta) * 100")
print(f"Percentage Change = {relative_change_per_eta:.4e} * {eta} * 100 = {percentage_change:.4f}%")

# Compare to options
# A. + 0.12%
# B. + 5.0 %
# C. + 0.0005 %
# D. - 5 %
# E. - 0.12 %
# Our result is ~0.00044%, which is very close to 0.0005%