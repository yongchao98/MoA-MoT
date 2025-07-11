import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power produced by Pu-239 after a specific
    uranium burnup percentage.
    """
    # Step 1: Define initial conditions and constants
    
    # Initial enrichment: 30% U-235, 70% U-238.
    # We can use these percentages as proportional initial numbers of atoms.
    N_U235_0 = 30.0
    N_U238_0 = 70.0
    N_U_0 = N_U235_0 + N_U238_0

    # Nuclear cross-sections in barns
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # Target burnup fraction
    burnup_fraction = 0.35
    target_U_consumed = burnup_fraction * N_U_0

    # Step 2: Define the equation to find the required neutron fluence (theta)
    def burnup_equation(theta):
        """
        This function represents the burnup condition. We need to find the root
        of this function, which is the fluence 'theta' that results in the
        target burnup.
        Equation: (Total Uranium Consumed) - (Target Uranium Consumed) = 0
        """
        # fsolve passes theta as an array
        t = theta[0]
        
        # Calculate the number of atoms of each uranium isotope consumed
        u235_consumed = N_U235_0 * (1 - np.exp(-sigma_a_U235 * t))
        u238_consumed = N_U238_0 * (1 - np.exp(-sigma_a_U238 * t))
        
        total_u_consumed = u235_consumed + u238_consumed
        
        return total_u_consumed - target_U_consumed

    # Solve for the neutron fluence 'theta' using a numerical solver
    # An initial guess of 0.03 is provided based on preliminary analysis.
    initial_guess_theta = [0.03]
    fluence_theta = fsolve(burnup_equation, initial_guess_theta)[0]

    # Step 3: Calculate the final nuclide concentrations at this fluence
    N_U235 = N_U235_0 * np.exp(-sigma_a_U235 * fluence_theta)
    
    # The concentration of Pu-239 is given by the solution to its rate equation
    # N_Pu239(t) = C * (exp(-sigma_a_U238*t) - exp(-sigma_a_Pu239*t))
    pu239_prefactor = (sigma_a_U238 * N_U238_0) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_u238 = np.exp(-sigma_a_U238 * fluence_theta)
    exp_term_pu239 = np.exp(-sigma_a_Pu239 * fluence_theta)
    N_Pu239 = pu239_prefactor * (exp_term_u238 - exp_term_pu239)

    # Step 4: Calculate the power contribution from each fissile isotope
    # Power is proportional to N * sigma_f
    power_U235 = N_U235 * sigma_f_U235
    power_Pu239 = N_Pu239 * sigma_f_Pu239

    total_power = power_U235 + power_Pu239
    fraction_from_Pu239 = power_Pu239 / total_power
    
    # Print the results, showing the numbers used in the final equation
    print("This script calculates the power fraction from Pu-239 at 35% total uranium burnup.\n")
    print(f"Calculated fluence (theta) corresponding to 35% burnup: {fluence_theta:.5f} [1/barns]")
    print("-" * 40)
    print("Final Proportional Concentrations:")
    print(f"N_U235  = {N_U235:.6f}")
    print(f"N_Pu239 = {N_Pu239:.6f}")
    print("-" * 40)
    print("Final Equation for Power Fraction from Pu-239:")
    print("Fraction = (N_Pu239 * sigma_f_Pu239) / (N_U235 * sigma_f_U235 + N_Pu239 * sigma_f_Pu239)")
    print(f"Fraction = ({N_Pu239:.6f} * {sigma_f_Pu239}) / ({N_U235:.6f} * {sigma_f_U235} + {N_Pu239:.6f} * {sigma_f_Pu239})")
    print(f"Fraction = ({power_Pu239:.6f}) / ({power_U235:.6f} + {power_Pu239:.6f})")
    print(f"Fraction = {fraction_from_Pu239:.6f}")

    return fraction_from_Pu239

# Run the calculation and store the final answer
final_answer = solve_power_fraction()
print(f"\nAfter 35% of the uranium has been burned, the fraction of power produced from Pu-239 is approximately {final_answer:.4f}, or {final_answer:.2%}.")
print(f"<<<{final_answer:.4f}>>>")
