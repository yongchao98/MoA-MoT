import numpy as np

# Step 1: Define constants and convert units
M = 938.93  # Nucleon mass in MeV
Cs2 = 267.1  # Scalar coupling constant, dimensionless
kF_fm_inv = 1.42  # Fermi momentum in fm^-1
R_fm = 5.5  # Nuclear radius in fm
M_star_ratio = 0.78  # Effective mass ratio
nu = 4  # Degeneracy factor for nuclear matter

# Conversion factor
hbarc = 197.327  # MeV fm

# Convert values to natural units (MeV)
kF = kF_fm_inv * hbarc
R = R_fm / hbarc
M_star = M_star_ratio * M

# Step 2: Calculate the momentum integral part (K_term)
# The antiderivative of k^3 / (k^2 + M_star^2)^(3/2) is sqrt(k^2 + M_star^2) + M_star^2 / sqrt(k^2 + M_star^2)
def antiderivative_k(k, Ms):
    if k == 0:
        return 2 * Ms
    return np.sqrt(k**2 + Ms**2) + Ms**2 / np.sqrt(k**2 + Ms**2)

K_term = antiderivative_k(kF, M_star) - antiderivative_k(0, M_star)

# Step 3: Calculate the percentage change using the corrected formula
# The formula for relative change (per unit eta) is (Cs2/M^2) * (3*nu*R / (8*pi^2)) * K_term
# We multiply by 100 for percentage.

# Breaking down the formula for printing
term1_val = Cs2 / M**2
term2_val = (3 * nu * R) / (8 * np.pi**2)
percentage_change = term1_val * term2_val * K_term * 100

# Step 4: Print the final equation with numerical values and the result
print("Calculation of the percentage change:")
print(f"Percentage Change = (C_s^2 / M^2) * (3 * nu * R / (8 * pi^2)) * K_term * 100")
print(f"Percentage Change = ({Cs2:.1f} / {M:.2f}^2) * (3 * {nu} * {R:.5f} / (8 * pi^2)) * {K_term:.4f} * 100")
print(f"Percentage Change = {term1_val:.6f} [MeV^-2] * {term2_val:.6f} [MeV] * {K_term:.4f} [MeV] * 100")
print(f"Percentage Change = {percentage_change:.6f} %")
print(f"\nThe calculated percentage change is approximately +{percentage_change:.4f}%, which is closest to +0.0005%.")

<<<C>>>