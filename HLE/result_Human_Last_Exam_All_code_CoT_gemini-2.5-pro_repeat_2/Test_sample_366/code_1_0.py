import numpy as np
from scipy.optimize import root_scalar

def solve_power_fraction():
    """
    Calculates the fraction of power produced by Plutonium-239
    after 35% of the total uranium fuel has been burned.
    """
    # --- Step 1: Define initial conditions and constants ---
    # Initial fuel composition (normalized)
    N_U235_0 = 0.30  # Initial fraction of U-235
    N_U238_0 = 0.70  # Initial fraction of U-238
    N_Pu239_0 = 0.0   # Initial fraction of Pu-239
    N_U_total_0 = N_U235_0 + N_U238_0

    # Absorption cross-sections (in barns)
    sigma_a_U235 = 591.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0

    # Fission cross-sections (in barns)
    sigma_f_U235 = 505.0
    sigma_f_Pu239 = 698.0

    # --- Step 2: Determine the neutron fluence (tau) ---
    # The condition is that 35% of total uranium is burned, so 65% remains.
    # We need to solve the following equation for tau:
    # N_U235_0*exp(-sigma_a_U235*tau) + N_U238_0*exp(-sigma_a_U238*tau) = 0.65 * N_U_total_0
    
    def burnup_equation(tau):
        u235_remaining = N_U235_0 * np.exp(-sigma_a_U235 * tau)
        u238_remaining = N_U238_0 * np.exp(-sigma_a_U238 * tau)
        return u235_remaining + u238_remaining - (0.65 * N_U_total_0)

    # Solve the equation numerically for tau. The root must be positive.
    # The solution τ will have units of barns^-1
    solution = root_scalar(burnup_equation, bracket=[0.0, 0.1], method='brentq')
    tau = solution.root
    
    # --- Step 3: Calculate final nuclide concentrations ---
    # Concentration of U-235 at fluence tau
    N_U235_f = N_U235_0 * np.exp(-sigma_a_U235 * tau)

    # Concentration of Pu-239 at fluence tau
    # This is the solution to the differential equation for Pu-239 production and destruction
    # dN_Pu/dτ = N_U238 * σ_a^U238 - N_Pu * σ_a^Pu
    term1 = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    term2 = np.exp(-sigma_a_U238 * tau) - np.exp(-sigma_a_Pu239 * tau)
    N_Pu239_f = term1 * term2

    # --- Step 4: Calculate the power fraction from Pu-239 ---
    # Power is proportional to N * sigma_f
    power_U235 = N_U235_f * sigma_f_U235
    power_Pu239 = N_Pu239_f * sigma_f_Pu239
    total_power = power_U235 + power_Pu239
    
    fraction_Pu239 = power_Pu239 / total_power

    # Print the final equation with the calculated values
    print("Calculation for the fraction of power from Plutonium-239:")
    print("Fraction = (N_Pu239 * σ_f^Pu239) / (N_U235 * σ_f^U235 + N_Pu239 * σ_f^Pu239)")
    print(f"Fraction = ({N_Pu239_f:.6f} * {sigma_f_Pu239}) / (({N_U235_f:.6e} * {sigma_f_U235}) + ({N_Pu239_f:.6f} * {sigma_f_Pu239}))")
    print("\nResult:")
    print(f"The fraction of power produced from Plutonium-239 is: {fraction_Pu239:.4f}")

    return fraction_Pu239

# Execute the function to get the answer
final_answer = solve_power_fraction()
print(f"<<<{final_answer:.4f}>>>")
