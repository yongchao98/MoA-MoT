import numpy as np
from scipy.optimize import root_scalar

# --- Step 1: Define initial conditions and constants ---

# Initial fraction of Uranium isotopes (based on 30% enrichment)
# We use these as the initial number of atoms in arbitrary units.
N_U235_0 = 0.30
N_U238_0 = 0.70
N_total_U_0 = N_U235_0 + N_U238_0

# Microscopic cross-sections in barns (1 barn = 1e-24 cm^2)
sigma_a_U235 = 591.0   # Absorption for U-235
sigma_f_U235 = 505.0   # Fission for U-235
sigma_a_U238 = 2.42    # Absorption for U-238
sigma_a_Pu239 = 973.0  # Absorption for Pu-239
sigma_f_Pu239 = 698.0  # Fission for Pu-239

print("Step 1: Initial Conditions and Constants")
print(f"Initial U-235 fraction (N_U235_0): {N_U235_0}")
print(f"Initial U-238 fraction (N_U238_0): {N_U238_0}")
print("-" * 30)

# --- Step 2: Solve for the Neutron Fluence (Φ) ---

# The problem states 35% of total uranium is burned, so 65% remains.
# N_U235(Φ) + N_U238(Φ) = 0.65 * N_total_U_0
# N_U235_0*exp(-σ_a_U235*Φ) + N_U238_0*exp(-σ_a_U238*Φ) = 0.65 * N_total_U_0
# We define a function whose root is the required fluence.
def find_fluence_equation(phi):
    """Equation for the remaining fraction of total uranium."""
    return (N_U235_0 * np.exp(-sigma_a_U235 * phi) +
            N_U238_0 * np.exp(-sigma_a_U238 * phi) - 0.65 * N_total_U_0)

# Numerically solve for the fluence (phi). We need a bracket [a, b] where
# the function values have opposite signs.
# f(0.01) is positive, f(0.05) is negative.
sol = root_scalar(find_fluence_equation, bracket=[0.01, 0.05], method='brentq')
fluence = sol.root

print("Step 2: Solve for Neutron Fluence (Φ)")
print(f"The neutron fluence (Φ) required to burn 35% of total uranium is: {fluence:.6f} barn^-1")
print("-" * 30)

# --- Step 3: Calculate final atom numbers ---

# Number of U-235 atoms remaining
N_U235 = N_U235_0 * np.exp(-sigma_a_U235 * fluence)

# Number of Pu-239 atoms produced and remaining.
# This is the solution to the Bateman equation for Pu-239 production and decay.
N_Pu239 = (N_U238_0 * sigma_a_U238 / (sigma_a_Pu239 - sigma_a_U238)) * \
          (np.exp(-sigma_a_U238 * fluence) - np.exp(-sigma_a_Pu239 * fluence))

print("Step 3: Calculate Final Atom Numbers at the given fluence")
print(f"Remaining U-235 atoms (N_U235): {N_U235:.4e}")
print(f"Remaining Pu-239 atoms (N_Pu239): {N_Pu239:.4e}")
print("-" * 30)

# --- Step 4: Calculate the power fraction ---

# Power is proportional to N * sigma_f
power_from_U235 = N_U235 * sigma_f_U235
power_from_Pu239 = N_Pu239 * sigma_f_Pu239
total_power = power_from_U235 + power_from_Pu239
fraction_from_Pu239 = power_from_Pu239 / total_power

print("Step 4: Calculate the Power Fraction from Plutonium-239")
print("\nThe final equation is: Fraction = Power_Pu239 / (Power_U235 + Power_Pu239)")
print("where Power is proportional to N * σ_f.\n")
print("Power from Pu-239 is proportional to:")
print(f"  N_Pu239 * σ_f^Pu239 = {N_Pu239:.4e} * {sigma_f_Pu239} = {power_from_Pu239:.4e}")
print("Power from U-235 is proportional to:")
print(f"  N_U235 * σ_f^U235  = {N_U235:.4e} * {sigma_f_U235} = {power_from_U235:.4e}")
print("\nFraction of power from Pu-239:")
print(f"  = {power_from_Pu239:.4e} / ({power_from_U235:.4e} + {power_from_Pu239:.4e})")
print(f"  = {fraction_from_Pu239:.6f}")

print(f"\nAfter 35% of the uranium has been burned, the fraction of power produced from plutonium-239 is {fraction_from_Pu239:.6f}.")
print("<<<" + f"{fraction_from_Pu239:.6f}" + ">>>")