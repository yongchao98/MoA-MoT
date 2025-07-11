import numpy as np
from scipy.optimize import fsolve

def solve_power_fraction():
    """
    Calculates the fraction of power produced from plutonium-239 in a nuclear reactor
    after a specific burnup percentage.
    """
    # Step 1: Define constants and initial conditions
    E = 0.30  # Initial U-235 enrichment
    B = 0.35  # Total uranium burnup fraction

    # Cross-sections in barns
    sigma_a_U235 = 591.0
    sigma_f_U235 = 505.0
    sigma_a_U238 = 2.42
    sigma_a_Pu239 = 973.0
    sigma_f_Pu239 = 698.0

    # Initial relative abundances (normalized to N_total_Uranium = 1)
    N_U235_0 = E
    N_U238_0 = 1 - E

    # Step 2: Define the function to find the neutron fluence (theta)
    def burnup_equation(theta):
        """
        Equation for burnup B as a function of fluence theta.
        We need to find the root of this function, i.e., where it equals zero.
        """
        # Using np.squeeze to ensure the output is a scalar if theta is an array
        theta = np.squeeze(theta)
        term_U235 = N_U235_0 * (1 - np.exp(-sigma_a_U235 * theta))
        term_U238 = N_U238_0 * (1 - np.exp(-sigma_a_U238 * theta))
        return term_U235 + term_U238 - B

    # Step 3: Solve for theta numerically
    # An initial guess is needed for the solver. A rough estimate suggests theta is small.
    initial_guess_theta = 0.01
    try:
        fluence_theta = fsolve(burnup_equation, initial_guess_theta)[0]
    except Exception as e:
        print(f"Error solving for fluence: {e}")
        return

    # Step 4: Calculate final isotope abundances using the solved fluence
    N_U235_final = N_U235_0 * np.exp(-sigma_a_U235 * fluence_theta)
    
    # The term (sigma_a_Pu239 - sigma_a_U238) is in the denominator
    denom = sigma_a_Pu239 - sigma_a_U238
    N_Pu239_final = (N_U238_0 * sigma_a_U238 / denom) * \
                    (np.exp(-sigma_a_U238 * fluence_theta) - np.exp(-sigma_a_Pu239 * fluence_theta))

    # Step 5: Calculate the power contributions and the final fraction
    power_from_U235 = N_U235_final * sigma_f_U235
    power_from_Pu239 = N_Pu239_final * sigma_f_Pu239
    
    total_power = power_from_U235 + power_from_Pu239
    
    # Avoid division by zero if total power is somehow zero
    if total_power == 0:
        fraction_from_Pu = 0
    else:
        fraction_from_Pu = power_from_Pu239 / total_power

    # Step 6: Output the results clearly
    print("The fraction of power from Plutonium-239 is calculated using the formula:")
    print("Fraction = (N_Pu239 * sigma_f_Pu239) / (N_U235 * sigma_f_U235 + N_Pu239 * sigma_f_Pu239)\n")
    
    print("Component values at 35% uranium burnup:")
    print(f"N_U235 (final relative abundance) = {N_U235_final:.6e}")
    print(f"sigma_f_U235 (barns) = {sigma_f_U235}")
    print(f"N_Pu239 (final relative abundance) = {N_Pu239_final:.6f}")
    print(f"sigma_f_Pu239 (barns) = {sigma_f_Pu239}\n")

    print("Plugging the numbers into the equation:")
    print(f"Fraction = ({N_Pu239_final:.6f} * {sigma_f_Pu239}) / (({N_U235_final:.6e} * {sigma_f_U235}) + ({N_Pu239_final:.6f} * {sigma_f_Pu239}))")
    print(f"Fraction = {power_from_Pu239:.6f} / ({power_from_U235:.6f} + {power_from_Pu239:.6f})")
    print(f"Fraction = {power_from_Pu239:.6f} / {total_power:.6f}\n")
    
    print(f"Final Answer: The fraction of power from Pu-239 is {fraction_from_Pu:.6f}")
    
    # Final answer in the required format
    global final_answer
    final_answer = fraction_from_Pu

if __name__ == '__main__':
    final_answer = 0.0
    solve_power_fraction()
    print(f"\n<<<_START_OF_ANSWER_>>>\n{final_answer}\n<<<END_OF_ANSWER>>>")
