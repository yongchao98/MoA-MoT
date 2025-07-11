import numpy as np
from scipy.optimize import root_scalar

def solve_power_fraction():
    """
    Calculates the fraction of power produced from plutonium-239
    after 35% of the initial uranium has been burned.
    """
    # --- Step 1: Define constants and initial conditions ---

    # Absorption cross-sections (in barns)
    sigma_a_U235 = 591.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0

    # Fission cross-sections (in barns)
    sigma_f_U235 = 505.0
    sigma_f_Pu239 = 698.0

    # Initial fuel composition (based on 30% enrichment, normalized to 100 atoms)
    N_U235_initial = 30.0
    N_U238_initial = 70.0
    N_U_initial = N_U235_initial + N_U238_initial
    N_Pu239_initial = 0.0
    
    # Burnup condition: 35% of uranium is burned, so 65% remains.
    N_U_final_target = 0.65 * N_U_initial

    # --- Step 2: Set up the equation for integrated flux ---

    # We need to find the integrated flux 'phi_t' such that the total number
    # of U atoms remaining matches the target.
    # The function below represents: N_U(t) - N_U_target = 0
    def uranium_balance(phi_t):
        n_u235_t = N_U235_initial * np.exp(-sigma_a_U235 * phi_t)
        n_u238_t = N_U238_initial * np.exp(-sigma_a_U238 * phi_t)
        return n_u235_t + n_u238_t - N_U_final_target

    # --- Step 3: Solve for the integrated flux (phi_t) ---
    
    # Use a numerical solver to find the root of the uranium_balance function.
    # The bracket [0.0, 0.1] provides a search range for the root.
    solution = root_scalar(uranium_balance, bracket=[0.0, 0.1])
    integrated_flux = solution.root

    # --- Step 4: Calculate final atom counts at the determined integrated flux ---
    
    # Number of U-235 atoms at time t
    N_U235_final = N_U235_initial * np.exp(-sigma_a_U235 * integrated_flux)

    # Number of Pu-239 atoms at time t, using the solution to the differential equation
    # N_Pu(t) = (sigma_a_U238 * N_U238_0) / (sigma_a_Pu239 - sigma_a_U238) * 
    #           (exp(-sigma_a_U238*t) - exp(-sigma_a_Pu239*t))
    common_term = (sigma_a_U238 * N_U238_initial) / (sigma_a_Pu239 - sigma_a_U238)
    exp_term_U238 = np.exp(-sigma_a_U238 * integrated_flux)
    exp_term_Pu239 = np.exp(-sigma_a_Pu239 * integrated_flux)
    N_Pu239_final = common_term * (exp_term_U238 - exp_term_Pu239)
    
    # --- Step 5: Calculate power components and the final fraction ---
    
    # Power is proportional to N * sigma_f
    power_U235 = N_U235_final * sigma_f_U235
    power_Pu239 = N_Pu239_final * sigma_f_Pu239
    
    total_power = power_U235 + power_Pu239
    
    # Final fraction of power from Plutonium-239
    fraction_from_Pu239 = power_Pu239 / total_power

    # Print the components of the final calculation
    print("Relative Power from U-235 Fission = {:.6f}".format(power_U235))
    print("Relative Power from Pu-239 Fission = {:.6f}".format(power_Pu239))
    print("\nFraction Equation:")
    print("Fraction = {:.6f} / ({:.6f} + {:.6f})".format(power_Pu239, power_U235, power_Pu239))

    print("\nFinal Answer:")
    print("The fraction of the power being produced from plutonium-239 is: {:.6f}".format(fraction_from_Pu239))
    
    # Returning the final answer in the required format
    return fraction_from_Pu239

# Run the calculation and print the result
final_fraction = solve_power_fraction()
print(f'<<<{final_fraction:.6f}>>>')
