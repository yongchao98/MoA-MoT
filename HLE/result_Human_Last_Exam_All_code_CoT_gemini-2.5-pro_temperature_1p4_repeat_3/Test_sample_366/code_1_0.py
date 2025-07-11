import numpy as np
from scipy.optimize import fsolve

# --- Step 1: Define Constants and Initial Conditions ---

# Cross-sections (in barns)
sigma_a_U235 = 591.0
sigma_f_U235 = 505.0
sigma_a_U238 = 2.42
sigma_a_Pu239 = 973.0
sigma_f_Pu239 = 698.0

# Initial relative abundances (as fractions).
# We can normalize the initial total number of Uranium atoms to 1.
N0_U235 = 0.30  # 30% enrichment
N0_U238 = 0.70  # 70% U-238

# Burnup fraction
burnup_fraction = 0.35

# --- Step 2: Solve for the Neutron Fluence (Phi) ---

# This function represents the burnup equation. We want to find the fluence (phi)
# that results in 35% of the initial uranium atoms being consumed.
# The equation is: (N0_U235 - N_final_U235) + (N0_U238 - N_final_U238) = 0.35 * (N0_U235 + N0_U238)
def burnup_equation(phi):
    # Number of U-235 atoms consumed
    burned_U235 = N0_U235 * (1 - np.exp(-sigma_a_U235 * phi))
    # Number of U-238 atoms consumed
    burned_U238 = N0_U238 * (1 - np.exp(-sigma_a_U238 * phi))
    # Total uranium atoms consumed
    total_burned = burned_U235 + burned_U238
    # Since N0_U235 + N0_U238 = 1, we want to find phi where total_burned = 0.35
    return total_burned - burnup_fraction

# Provide an initial guess for the solver and find the root.
initial_guess_phi = 0.01
fluence_phi = fsolve(burnup_equation, initial_guess_phi)[0]


# --- Step 3: Calculate Final Nuclide Concentrations ---

# Calculate the final relative abundance of U-235
N_final_U235 = N0_U235 * np.exp(-sigma_a_U235 * fluence_phi)

# Calculate the final relative abundance of Pu-239 using the Bateman equation for the U-238 -> Pu-239 chain.
# N_Pu239(t) = N0_U238 * (lambda_U238 / (lambda_Pu239 - lambda_U238)) * (exp(-lambda_U238*t) - exp(-lambda_Pu239*t))
# where lambda = sigma_a * flux. The flux cancels out.
pre_factor = N0_U238 * (sigma_a_U238 / (sigma_a_Pu239 - sigma_a_U238))
exp_term_diff = np.exp(-sigma_a_U238 * fluence_phi) - np.exp(-sigma_a_Pu239 * fluence_phi)
N_final_Pu239 = pre_factor * exp_term_diff


# --- Step 4: Calculate and Print the Power Fraction ---

print("This script calculates the fraction of power from Plutonium-239 after 35% uranium burnup.\n")

print("--- Calculation Steps ---")
print(f"1. Solved for neutron fluence (Φ) to achieve 35% burnup: {fluence_phi:.6f} [1/barns]")

print(f"\n2. Calculated final relative abundances (normalized to initial total Uranium atoms):")
print(f"   - Final U-235: {N_final_U235:.6f}")
print(f"   - Final Pu-239: {N_final_Pu239:.6f}")

# Power is proportional to N * sigma_f. The proportionality constants cancel in the fraction.
power_from_U235 = N_final_U235 * sigma_f_U235
power_from_Pu239 = N_final_Pu239 * sigma_f_Pu239
total_power = power_from_U235 + power_from_Pu239
fraction_power_Pu239 = power_from_Pu239 / total_power

print("\n3. Calculated the components of the power fraction at the final state:")
print("   (Power is proportional to N_final * σ_f)")
print(f"   - Power term from Pu-239 = N_final_Pu239 * σ_f_Pu239 = {N_final_Pu239:.6f} * {sigma_f_Pu239} = {power_from_Pu239:.4f}")
print(f"   - Power term from U-235 = N_final_U235 * σ_f_U235 = {N_final_U235:.6f} * {sigma_f_U235} = {power_from_U235:.4f}")

print("\n4. Final power fraction calculation:")
print(f"   Fraction = (Power from Pu-239) / (Power from Pu-239 + Power from U-235)")
print(f"   Fraction = {power_from_Pu239:.4f} / ({power_from_Pu239:.4f} + {power_from_U235:.4f})")
print(f"   Fraction = {power_from_Pu239:.4f} / {total_power:.4f}")

print("\n---------------------------------------------------------------------------------")
print(f"The fraction of power being produced from plutonium-239 is: {fraction_power_Pu239:.4f}")
print("---------------------------------------------------------------------------------")

print(f"\n<<< {fraction_power_Pu239:.4f} >>>")