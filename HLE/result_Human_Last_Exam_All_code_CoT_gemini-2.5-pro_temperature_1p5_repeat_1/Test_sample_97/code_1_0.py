import math

def calculate_solubility():
    """
    Calculates the solubility of Al(OH)3 in pure water, considering both
    simple dissolution and complex ion formation.
    """
    # Known constants
    K_sp = 5.3e-27  # Solubility product of Al(OH)3
    K_f = 1.1e31   # Formation constant of Al(OH)4-
    K_w = 1.0e-14  # Ion product of water

    # The governing equation for the system is derived from the charge balance:
    # 3[Al^3+] + [H+] = [OH-] + [Al(OH)4-]
    # Substituting equilibrium expressions leads to a polynomial in [OH-]:
    # (1 + K_sp*K_f)*[OH-]^4 - K_w*[OH-]^2 - 3*K_sp = 0
    # This is a quadratic equation in y = [OH-]^2, of the form ay^2 + by + c = 0.

    # Coefficients of the quadratic equation for y = [OH-]^2
    a = 1 + K_sp * K_f
    b = -K_w
    c = -3 * K_sp

    # Solve the quadratic equation for y. We take the positive root because y represents
    # a physical quantity ([OH-]^2) which must be positive.
    discriminant = b**2 - 4 * a * c
    y = (-b + math.sqrt(discriminant)) / (2 * a)

    # [OH-] is the square root of y
    oh_concentration = math.sqrt(y)

    # Calculate the concentrations of the individual aluminum species at equilibrium
    al3_concentration = K_sp / (oh_concentration**3)
    al_oh4_concentration = K_sp * K_f * oh_concentration

    # The total solubility is the sum of the concentrations of all dissolved aluminum species
    solubility = al3_concentration + al_oh4_concentration

    print("Step 1: Calculate the equilibrium hydroxide concentration [OH-].")
    print(f"Solving the governing polynomial equation yields [OH-] = {oh_concentration:.2e} mol L^-1.\n")

    print("Step 2: Calculate the concentrations of the dissolved aluminum species.")
    print(f"[Al^3+] = K_sp / [OH-]^3 = {al3_concentration:.2e} mol L^-1")
    print(f"[Al(OH)4-] = K_sp * K_f * [OH-] = {al_oh4_concentration:.2e} mol L^-1\n")

    print("Step 3: Calculate the total solubility (S).")
    print("S = [Al^3+] + [Al(OH)4-]")
    # The final equation with the calculated numbers
    print(f"S = {al3_concentration:.2e} mol L^-1 + {al_oh4_concentration:.2e} mol L^-1")
    
    # Final answer, rounded to two significant figures consistent with the input constants
    print(f"S = {solubility:.2e} mol L^-1")

calculate_solubility()