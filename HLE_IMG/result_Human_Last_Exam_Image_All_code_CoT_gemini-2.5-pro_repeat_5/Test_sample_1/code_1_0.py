import math

# --- Step 1: Define constants and energies from the problem ---
# Constants
k_B = 1.380649e-23  # Boltzmann constant in J/K
h = 6.62607015e-34   # Planck constant in J*s
R = 8.31446         # Ideal gas constant in J/(mol*K)
T = 298             # Temperature in K

# Gibbs free energies from the diagram in kJ/mol
G_R = -55
G_TS_add = 10
G_P5 = -22
G_TS_eli = 18

# --- Step 2: Identify the Rate-Determining Step (RDS) ---
# Calculate the activation barrier for each chemical transformation step
barrier_1 = G_TS_add - G_R  # Barrier for R -> P5
barrier_2 = G_TS_eli - G_P5  # Barrier for P5 -> P

print("--- Analysis of Rate-Determining Step ---")
print(f"Activation barrier for water addition (R -> P5): {G_TS_add} - ({G_R}) = {barrier_1} kJ/mol")
print(f"Activation barrier for fluoride elimination (P5 -> P): {G_TS_eli} - ({G_P5}) = {barrier_2} kJ/mol")

# The RDS is the one with the highest activation barrier
if barrier_1 > barrier_2:
    rds_description = "the water addition step (R -> P5)"
    delta_G_ddagger_kJ = barrier_1
else:
    rds_description = "the fluoride elimination step (P5 -> P)"
    delta_G_ddagger_kJ = barrier_2

print(f"\nThe rate-determining step is {rds_description} with an activation barrier of {delta_G_ddagger_kJ} kJ/mol.\n")

# --- Step 3: Calculate the Reaction Rate Constant ---
# Convert activation energy to J/mol
delta_G_ddagger_J = delta_G_ddagger_kJ * 1000

# Calculate the pre-exponential factor (k_B * T / h)
pre_exponential_factor = (k_B * T) / h

# Calculate the exponential term exp(-ΔG‡ / RT)
exponent = -delta_G_ddagger_J / (R * T)
exponential_factor = math.exp(exponent)

# Calculate the rate constant in s^-1
k_seconds = pre_exponential_factor * exponential_factor

# Convert the rate constant to hours^-1
k_hours = k_seconds * 3600

# Round to two significant figures
k_hours_rounded = float(f"{k_hours:.2g}")

print("--- Calculation of the Rate Constant (k) ---")
print("Using the Eyring equation: k = (k_B * T / h) * exp(-ΔG‡ / (R * T))")
print(f"k = ({k_B:.4g} J/K * {T} K / {h:.4g} J*s) * exp(-{delta_G_ddagger_J} J/mol / ({R:.4g} J/(mol*K) * {T} K))")
print(f"k = ({pre_exponential_factor:.3e} s^-1) * exp({exponent:.3f})")
print(f"k = ({pre_exponential_factor:.3e} s^-1) * ({exponential_factor:.3e})")
print(f"k = {k_seconds:.2f} s^-1")
print(f"k = {k_seconds:.2f} s^-1 * 3600 s/hr = {k_hours:.2f} hr^-1")
print(f"\nThe reaction rate constant to two significant figures is {k_hours_rounded:.2g} hours^-1.")

print(f"\n<<<9.1e+04>>>")