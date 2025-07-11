import math

def calculate_solubility():
    """
    Calculates the solubility of Al(OH)3 in pure water, considering both dissolution
    and complex ion formation.
    """
    # Known constants
    K_sp = 5.3e-27
    K_f = 1.1e-31
    K_w = 1.0e-14
    
    # K for the reaction Al(OH)3(s) + OH-(aq) <=> Al(OH)4-(aq)
    K_c = K_sp * K_f

    # The charge balance equation leads to a quadratic equation for x = [OH-]^2:
    # (1 + K_sp*K_f)*[OH-]^4 - K_w*[OH-]^2 - 3*K_sp = 0
    # Let x = [OH-]^2, then we have Ax^2 + Bx + C = 0
    
    A = 1 + K_c
    B = -K_w
    C = -3 * K_sp

    # Solve the quadratic equation for x = [OH-]^2 using the quadratic formula
    # x = (-B + sqrt(B^2 - 4AC)) / 2A (we take the positive root)
    discriminant = math.sqrt(B**2 - 4*A*C)
    x = (-B + discriminant) / (2*A)

    # Calculate concentrations
    oh_conc = math.sqrt(x)
    al3_conc = K_sp / (oh_conc**3)
    aloh4_conc = K_c * oh_conc
    
    # Total solubility S = [Al^3+] + [Al(OH)4^-]
    solubility = al3_conc + aloh4_conc

    # Print the results
    print("This script calculates the molar solubility of Al(OH)3 in pure water.")
    print("The calculation is based on charge balance and simultaneous equilibria.\n")
    print(f"The equilibrium hydroxide concentration [OH-] is calculated to be: {oh_conc:.4e} mol L^-1\n")
    print("The final solubility (S) is the sum of the dissolved aluminum species:")
    print("S = [Al^3+] + [Al(OH)4^-]\n")
    
    print("The concentration of Al^3+ is:")
    print(f"[Al^3+] = K_sp / [OH-]^3 = {K_sp:.2e} / ({oh_conc:.4e})^3 = {al3_conc:.4e} mol L^-1")
    
    print("\nThe concentration of Al(OH)4^- is:")
    print(f"[Al(OH)4^-] = (K_sp * K_f) * [OH-] = ({K_c:.2e}) * {oh_conc:.4e} = {aloh4_conc:.4e} mol L^-1")

    print("\n----------------------------------------------------------")
    print("Total Molar Solubility:")
    print(f"S = {al3_conc:.4e} + {aloh4_conc:.4e} = {solubility:.4e} mol L^-1")
    print("----------------------------------------------------------")
    
    # Return the final answer in the required format
    print(f"\n<<<{solubility:.4e}>>>")

calculate_solubility()