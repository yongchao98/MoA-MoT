import numpy as np
from scipy.optimize import fsolve

# Step 1: Define physical constants and initial conditions
s_a_U235 = 591.0   # Absorption cross-section for U-235 (barns)
s_f_U235 = 505.0   # Fission cross-section for U-235 (barns)
s_a_U238 = 2.42    # Absorption cross-section for U-238 (barns)
s_a_Pu239 = 973.0  # Absorption cross-section for Pu-239 (barns)
s_f_Pu239 = 698.0  # Fission cross-section for Pu-239 (barns)

# Initial fuel composition (fractional)
N0_U235_frac = 0.30
N0_U238_frac = 0.70

# Condition: 35% of total uranium is burned, so 65% remains.
remaining_U_frac = 0.65

# Step 2: Numerically solve for the neutron fluence (Phi)
# We need to find the root of the function f(phi) = 0, where:
# f(phi) = (U-235 remaining) + (U-238 remaining) - (target remaining fraction)
def find_fluence_eq(phi):
    """Equation to find the fluence that results in the target uranium remaining fraction."""
    return (N0_U235_frac * np.exp(-s_a_U235 * phi) +
            N0_U238_frac * np.exp(-s_a_U238 * phi) -
            remaining_U_frac)

# An initial guess for the root-finding algorithm.
# Based on preliminary analysis, the fluence should be around 0.03.
initial_guess_phi = 0.03
# Use fsolve to find the fluence. We take the first element of the result array.
fluence = fsolve(find_fluence_eq, initial_guess_phi)[0]

# Step 3: Calculate the number densities of U-235 and Pu-239 at this fluence.
# These are fractions relative to the initial total number of uranium atoms.
N_U235 = N0_U235_frac * np.exp(-s_a_U235 * fluence)

# Formula for Pu-239 concentration
term1 = (N0_U238_frac * s_a_U238) / (s_a_Pu239 - s_a_U238)
term2 = np.exp(-s_a_U238 * fluence) - np.exp(-s_a_Pu239 * fluence)
N_Pu239 = term1 * term2

# Step 4: Calculate the power fraction from Pu-239 fission
# Power is proportional to the macroscopic fission rate (N * sigma_f)
power_U235 = N_U235 * s_f_U235
power_Pu239 = N_Pu239 * s_f_Pu239
total_power = power_U235 + power_Pu239
power_fraction_Pu239 = power_Pu239 / total_power

# Print the intermediate values and the final result
print(f"Calculation Steps for Power Fraction from Plutonium-239:\n")
print(f"1. Calculated Neutron Fluence (Φ) for 35% Uranium Burnup: {fluence:.6f} [1/barn]")
print(f"2. Relative atom densities at this fluence:")
print(f"   N_U235 = {N_U235:.6e} (relative to initial total U atoms)")
print(f"   N_Pu239 = {N_Pu239:.6f} (relative to initial total U atoms)\n")

print("3. Power fraction is calculated as: (N_Pu239 * σ_f_Pu239) / (N_U235 * σ_f_U235 + N_Pu239 * σ_f_Pu239)")
print(f"   Fraction = ({N_Pu239:.6f} * {s_f_Pu239}) / (({N_U235:.6e} * {s_f_U235}) + ({N_Pu239:.6f} * {s_f_Pu239}))")
print(f"   Fraction = {power_Pu239:.6f} / ({power_U235:.6e} + {power_Pu239:.6f})")
print(f"   Fraction = {power_Pu239:.6f} / {total_power:.6f}\n")
print(f"Final Answer: The fraction of power being produced from plutonium-239 is {power_fraction_Pu239:.4f}")
print(f"<<<{power_fraction_Pu239:.4f}>>>")