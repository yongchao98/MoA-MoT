import numpy as np
from scipy.optimize import root_scalar

# Plan Step 1: Define initial conditions and constants from the problem.
# We normalize the initial total number of Uranium atoms to 1.
N0_U235 = 0.30  # Initial fraction of U-235
N0_U238 = 0.70  # Initial fraction of U-238
N0_Total_U = N0_U235 + N0_U238

# Absorption (a) and fission (f) cross-sections (sigma) in barns
sa_U235 = 591.0
sf_U235 = 505.0
sa_U238 = 2.42
sa_Pu239 = 973.0
sf_Pu239 = 698.0

# Plan Step 2: Define the burnup equation and solve for fluence (tau).
# The equation describes the remaining fraction of total uranium.
# We need to find tau such that N_U235(tau) + N_U238(tau) = 0.65 * N0_Total_U
def uranium_burnup_equation(tau):
    """
    Calculates the difference between the remaining uranium fraction and the target.
    We are looking for the root of this function.
    """
    n_u235_t = N0_U235 * np.exp(-sa_U235 * tau)
    n_u238_t = N0_U238 * np.exp(-sa_U238 * tau)
    # Target condition: 35% of uranium burned means 65% remains.
    return (n_u235_t + n_u238_t) - (0.65 * N0_Total_U)

# Solve for tau, the neutron fluence in units of barns^-1
# We'll search for a root in the interval [0.01, 0.1] based on preliminary analysis.
sol = root_scalar(uranium_burnup_equation, bracket=[0.01, 0.1], method='brentq')
tau = sol.root

print("--- Step 1 & 2: Determine Neutron Fluence (τ) for 35% Uranium Burnup ---")
print(f"Initial U-235 fraction (N0_U235): {N0_U235}")
print(f"Initial U-238 fraction (N0_U238): {N0_U238}")
print(f"Solved neutron fluence (τ) for 35% burnup: {tau:.5f} barns^-1")
print("-" * 60)

# Plan Step 3: Calculate the number densities of each isotope at the calculated fluence.
# Equation for N_U235 at fluence tau
N_U235 = N0_U235 * np.exp(-sa_U235 * tau)
# Equation for N_U238 at fluence tau
N_U238 = N0_U238 * np.exp(-sa_U238 * tau)
# Equation for N_Pu239 at fluence tau
# N_Pu239(t) = N0_U238*sa_U238/(sa_Pu239 - sa_U238) * (exp(-sa_U238*t) - exp(-sa_Pu239*t))
N_Pu239 = (N0_U238 * sa_U238) / (sa_Pu239 - sa_U238) * (np.exp(-sa_U238 * tau) - np.exp(-sa_Pu239 * tau))

print("--- Step 3: Calculate Isotope Amounts at Target Fluence ---")
print(f"Relative amount of U-235 (N_U235): {N_U235:.4e}")
print(f"Relative amount of U-238 (N_U238): {N_U238:.4e}")
print(f"Relative amount of Pu-239 (N_Pu239): {N_Pu239:.4e}")
# Verification of burnup
print(f"Verification: Total remaining U / Initial U = {(N_U235 + N_U238) / N0_Total_U:.3f}")
print("-" * 60)

# Plan Step 4: Calculate the power fraction from Pu-239.
# Power is proportional to the fission rate (R = N * sigma_f)
# Note: The neutron flux term (phi) cancels out in the fraction.
Power_U235 = N_U235 * sf_U235
Power_Pu239 = N_Pu239 * sf_Pu239
Total_Power = Power_U235 + Power_Pu239

if Total_Power > 0:
    fraction_Pu_power = Power_Pu239 / Total_Power
else:
    fraction_Pu_power = 0

print("--- Step 4: Calculate Power Fraction from Plutonium-239 ---")
print("Power is proportional to N * σ_f.")
print(f"Power from U-235 is proportional to {N_U235:.4e} * {sf_U235} = {Power_U235:.4e}")
print(f"Power from Pu-239 is proportional to {N_Pu239:.4e} * {sf_Pu239} = {Power_Pu239:.4e}")
print("\nThe fraction of power from Plutonium is:")
print(f"Fraction = Power_Pu239 / (Power_U235 + Power_Pu239)")
print(f"Fraction = {Power_Pu239:.4f} / ({Power_U235:.4f} + {Power_Pu239:.4f})")
print(f"Fraction = {Power_Pu239:.4f} / {Total_Power:.4f}")
print("\nFinal Answer:")
print(f"The fraction of power produced from Plutonium-239 is: {fraction_Pu_power:.4f}")
print(f"<<<{fraction_Pu_power:.4f}>>>")