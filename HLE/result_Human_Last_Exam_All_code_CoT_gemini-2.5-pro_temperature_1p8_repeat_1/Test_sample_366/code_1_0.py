import numpy as np
from scipy.optimize import root_scalar

# Step 1: Define System Parameters
# Cross-sections in barns
sigma_a_U235 = 591.0
sigma_f_U235 = 505.0
sigma_a_U238 = 2.42
sigma_a_Pu239 = 973.0
sigma_f_Pu239 = 698.0

# Initial fractions (atom fractions)
f_U235_initial = 0.30
f_U238_initial = 0.70

# Total uranium burnup fraction
burnup_fraction = 0.35

# We use relative initial number densities, assuming total initial Uranium atoms N_U_initial = 1.0
N_U235_initial = f_U235_initial
N_U238_initial = f_U238_initial

# Step 2: Define and Solve the Burnup Equation for Neutron Fluence (theta)
# The problem states 35% of total Uranium is burned. This means 65% remains.
# The total amount of uranium remaining is the sum of U-235 and U-238 remaining.
# We must find the neutron fluence (theta) that satisfies this condition.
def burnup_equation(theta):
    """
    This function represents the state of the fuel. We need to find the root, theta,
    where the final amount of uranium equals the initial amount minus the burned amount.
    The equation is: (U235_final + U238_final) - (1 - burnup_fraction) = 0
    """
    final_U235 = N_U235_initial * np.exp(-sigma_a_U235 * theta)
    final_U238 = N_U238_initial * np.exp(-sigma_a_U238 * theta)
    return final_U235 + final_U238 - (1.0 - burnup_fraction)

# Step 3: Solve for theta using a numerical root finder
# We search for a root in a plausible bracket [0.0, 0.1].
sol = root_scalar(burnup_equation, bracket=[0.0, 0.1], method='brentq')
theta = sol.root

# Step 4: Calculate Final Isotope Concentrations
# With the fluence (theta), find the final relative number of atoms for each isotope.
N_U235_final = N_U235_initial * np.exp(-sigma_a_U235 * theta)

# The concentration of Pu-239 is given by the Bateman equation solution.
# It is produced from U-238 capture and destroyed by its own neutron absorption.
N_Pu239_final = (N_U238_initial * sigma_a_U238 / (sigma_a_Pu239 - sigma_a_U238)) * \
                (np.exp(-sigma_a_U238 * theta) - np.exp(-sigma_a_Pu239 * theta))

# Step 5: Calculate Power Fraction from Plutonium-239
# Power is proportional to the fission rate R_f = N * sigma_f.
# The neutron flux cancels out when taking the ratio.

# Power term from U-235 fission
power_U235 = N_U235_final * sigma_f_U235

# Power term from Pu-239 fission
power_Pu239 = N_Pu239_final * sigma_f_Pu239

# Total fission power term
total_power = power_U235 + power_Pu239

# Fraction of power from Pu-239
fraction_power_Pu239 = power_Pu239 / total_power

# Step 6: Output the calculation and final result
print("The fraction of power produced by plutonium is calculated as the ratio of the fission rate of plutonium to the total fission rate.")
print("Fraction = (Power from Pu-239) / (Total Power)")
print("         = (N_Pu239 * sigma_f_Pu239) / ( (N_U235 * sigma_f_U235) + (N_Pu239 * sigma_f_Pu239) )")
print("\nAfter solving for the neutron fluence required for 35% uranium burnup, we find the final concentrations and power terms.")
print(f"\nFinal relative concentration of U-235 (N_U235): {N_U235_final:.5e}")
print(f"Final relative concentration of Pu-239 (N_Pu239): {N_Pu239_final:.5f}")
print(f"\nFission cross-section for U-235 (sigma_f_U235): {sigma_f_U235} barns")
print(f"Fission cross-section for Pu-239 (sigma_f_Pu239): {sigma_f_Pu239} barns")
print("\nThe final equation with the calculated values is:")
print(f"Fraction = ({power_Pu239:.5f}) / ( ({power_U235:.5e}) + ({power_Pu239:.5f}) )")
print(f"Fraction = {power_Pu239:.5f} / ({total_power:.5f})")

print(f"\nThe final fraction of power from Pu-239 is: {fraction_power_Pu239:.5f}")

<<<0.99999>>>