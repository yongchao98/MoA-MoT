import math

def calculate_solubility():
    """
    Calculates the solubility of Al(OH)3 in pure water,
    considering the formation of the Al(OH)4- complex ion.
    """
    # Known constants
    K_sp = 5.3e-27
    K_f = 1.1e31
    K_w = 1.0e-14

    # The charge balance equation leads to a quadratic equation for x = [OH-]^2:
    # a*x^2 + b*x + c = 0
    # where x = [OH-]^2
    a = 1 + K_sp * K_f
    b = -K_w
    c = -3 * K_sp

    # Solve the quadratic equation for x = [OH-]^2
    # x = (-b + sqrt(b^2 - 4ac)) / (2a)
    # We take the positive root since concentration squared must be positive.
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("Error: No real solution for [OH-]. Check constants.")
        return

    oh_conc_squared = (-b + math.sqrt(discriminant)) / (2 * a)
    oh_conc = math.sqrt(oh_conc_squared)

    # Calculate the concentrations of the aluminum species
    al3_conc = K_sp / (oh_conc**3)
    al_oh4_conc = K_sp * K_f * oh_conc

    # The total solubility S is the sum of the concentrations of dissolved aluminum species
    solubility = al3_conc + al_oh4_conc

    # Print the breakdown of the final solubility calculation
    print(f"The calculation for solubility (S) is the sum of the concentrations of the aluminum species:")
    print(f"S = [Al^3+] + [Al(OH)4^-]")
    print(f"S = {al3_conc:.3e} mol L^-1 + {al_oh4_conc:.3e} mol L^-1")
    
    # Print the final answer
    print(f"\nThe total solubility of Al(OH)3 in pure water is {solubility:.2e} mol L^-1.")

calculate_solubility()
<<<1.8e-3>>>