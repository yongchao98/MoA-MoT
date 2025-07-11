import math
from scipy.optimize import fsolve

def solve_reactor_burnup():
    """
    Calculates the fraction of power produced from Pu-239 after a specific
    uranium burnup percentage.
    """
    # Step 1: Define constants and initial conditions
    # Cross-sections in barns (1 barn = 1e-24 cm^2)
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # Initial enrichment and burnup fraction
    initial_enrichment = 0.30
    burnup_fraction = 0.35

    # Set a basis for the initial number of atoms
    N_U_0 = 100.0
    N_U235_0 = initial_enrichment * N_U_0
    N_U238_0 = (1 - initial_enrichment) * N_U_0

    # Step 2: Define the equation to find the neutron fluence (theta)
    # We need to find theta such that the remaining uranium is (1 - burnup_fraction)
    # of the initial amount.
    def find_fluence_equation(theta, N_U235_0, N_U238_0, N_U_0, burnup_fraction, sigma_a_U235, sigma_a_U238):
        remaining_U235 = N_U235_0 * math.exp(-sigma_a_U235 * theta)
        remaining_U238 = N_U238_0 * math.exp(-sigma_a_U238 * theta)
        target_N_U = (1 - burnup_fraction) * N_U_0
        return remaining_U235 + remaining_U238 - target_N_U

    # Step 3: Solve for the fluence theta numerically
    # An initial guess for theta is needed for the solver.
    initial_guess_theta = 0.001
    args = (N_U235_0, N_U238_0, N_U_0, burnup_fraction, sigma_a_U235, sigma_a_U238)
    fluence_theta = fsolve(find_fluence_equation, initial_guess_theta, args=args)[0]

    # Step 4: Calculate the final number of atoms for each nuclide
    N_U235_final = N_U235_0 * math.exp(-sigma_a_U235 * fluence_theta)
    
    # Equation for Pu-239 production and destruction
    term_pu_1 = (N_U238_0 * sigma_a_U238) / (sigma_a_Pu239 - sigma_a_U238)
    term_pu_2 = math.exp(-sigma_a_U238 * fluence_theta) - math.exp(-sigma_a_Pu239 * fluence_theta)
    N_Pu239_final = term_pu_1 * term_pu_2

    # Step 5: Calculate the power fraction from Pu-239
    # Power is proportional to N * sigma_f
    power_U235 = N_U235_final * sigma_f_U235
    power_Pu239 = N_Pu239_final * sigma_f_Pu239
    
    total_power = power_U235 + power_Pu239
    power_fraction_Pu239 = power_Pu239 / total_power

    # Print the results and the final calculation steps
    print("--- Calculation Details ---")
    print(f"Initial relative number of U-235 atoms: {N_U235_0}")
    print(f"Initial relative number of U-238 atoms: {N_U238_0}")
    print(f"Solved neutron fluence (theta): {fluence_theta:.6f} [1/barns]")
    print("\n--- Final Nuclide Concentrations (Relative) ---")
    print(f"Final number of U-235 atoms: {N_U235_final:.4f}")
    print(f"Final number of Pu-239 atoms: {N_Pu239_final:.4f}")
    
    print("\n--- Final Power Fraction Calculation ---")
    print(f"The fraction of power from Plutonium-239 is given by:")
    print(f"Fraction = (N_Pu239 * σ_f_Pu239) / (N_U235 * σ_f_U235 + N_Pu239 * σ_f_Pu239)\n")
    
    print(f"Power from Pu-239 is proportional to: {N_Pu239_final:.4f} * {sigma_f_Pu239} = {power_Pu239:.2f}")
    print(f"Power from U-235 is proportional to: {N_U235_final:.4f} * {sigma_f_U235} = {power_U235:.2f}")
    
    print(f"\nFraction = {power_Pu239:.2f} / ({power_U235:.2f} + {power_Pu239:.2f})")
    print(f"Fraction = {power_Pu239:.2f} / {total_power:.2f}")
    
    print(f"\nFinal Answer: The fraction of power from Pu-239 is {power_fraction_Pu239:.4f}")
    
    return power_fraction_Pu239

# Execute the function and store the final answer
final_answer = solve_reactor_burnup()
print(f"\n<<<{final_answer:.4f}>>>")
