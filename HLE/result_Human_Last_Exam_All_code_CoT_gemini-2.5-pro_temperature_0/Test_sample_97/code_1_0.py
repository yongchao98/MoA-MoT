import math

def solve_solubility():
    """
    Calculates the solubility of Al(OH)3 in pure water, considering complex formation.
    """
    # Given constants
    K_sp = 5.3e-27
    K_f = 1.1e31
    K_w = 1.0e-14

    # Calculate the equilibrium constant for the complex formation from the solid
    # Al(OH)3(s) + OH-(aq) <=> Al(OH)4-(aq)
    # K = K_sp * K_f
    K = K_sp * K_f

    print("Step 1: Define constants and set up the main equation.")
    print(f"K_sp = {K_sp:.2e}")
    print(f"K_f = {K_f:.2e}")
    print(f"K_w = {K_w:.2e}")
    print(f"The combined constant K = K_sp * K_f = {K:.2e}")
    print("\nThe charge balance equation 3[Al^3+] + [H+] = [OH^-] + [Al(OH)4^-] leads to:")
    print("(K+1)*[OH-]^4 - K_w*[OH-]^2 - 3*K_sp = 0")
    print("Let y = [OH-]^2, the equation becomes a quadratic: ay^2 + by + c = 0\n")

    # Set up the quadratic equation for y = [OH-]^2
    # a*y^2 + b*y + c = 0
    a = K + 1
    b = -K_w
    c = -3 * K_sp

    print("Step 2: Solve the quadratic equation for y = [OH-]^2.")
    print(f"a = {a:.2e}")
    print(f"b = {b:.2e}")
    print(f"c = {c:.2e}\n")

    # Solve using the quadratic formula: y = [-b ± sqrt(b^2 - 4ac)] / 2a
    discriminant = math.sqrt(b**2 - 4*a*c)
    # Since y = [OH-]^2 must be positive, we take the positive root
    y = (-b + discriminant) / (2*a)

    # Calculate [OH-]
    hydroxide_conc = math.sqrt(y)
    print(f"Step 3: Calculate the equilibrium hydroxide concentration [OH-].")
    print(f"[OH-] = {hydroxide_conc:.3e} mol L^-1\n")

    # Calculate the concentrations of the aluminum species
    al_ion_conc = K_sp / (hydroxide_conc**3)
    complex_ion_conc = K * hydroxide_conc
    
    print("Step 4: Calculate the concentrations of dissolved aluminum species.")
    print(f"[Al^3+] = K_sp / [OH-]^3 = {K_sp:.2e} / ({hydroxide_conc:.3e})^3 = {al_ion_conc:.3e} mol L^-1")
    print(f"[Al(OH)4^-] = K * [OH-] = {K:.2e} * {hydroxide_conc:.3e} = {complex_ion_conc:.3e} mol L^-1\n")

    # Calculate the total solubility
    solubility = al_ion_conc + complex_ion_conc
    
    print("Step 5: Calculate the total solubility (S).")
    print(f"S = [Al^3+] + [Al(OH)4^-]")
    print(f"S = {al_ion_conc:.3e} + {complex_ion_conc:.3e}")
    print(f"The total solubility of Al(OH)3 is {solubility:.3e} mol L^-1.")
    
    # Return the final answer in the required format
    print(f"\n<<<{solubility:.2e}>>>")

solve_solubility()