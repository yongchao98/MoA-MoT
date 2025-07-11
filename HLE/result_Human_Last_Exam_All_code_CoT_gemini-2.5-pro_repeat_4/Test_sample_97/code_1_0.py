import math

def calculate_solubility():
    """
    Calculates the solubility of Al(OH)3 in pure water,
    considering the formation of the Al(OH)4- complex ion.
    """
    # --- Step 1: Define constants ---
    K_sp = 5.3e-27
    K_f = 1.1e31
    K_w = 1.0e-14

    print("--- Calculation of Al(OH)3 Solubility ---")
    print("\nStep 1: Define the equilibrium constants.")
    print(f"K_sp for Al(OH)3 <=> Al^3+ + 3OH^- is {K_sp:.2e}")
    print(f"K_f for Al^3+ + 4OH^- <=> Al(OH)4^- is {K_f:.2e}")
    print(f"K_w for H2O <=> H+ + OH- is {K_w:.2e}")

    # --- Step 2: Calculate the combined equilibrium constant K_c ---
    K_c = K_sp * K_f
    print("\nStep 2: Combine equilibria to relate Al(OH)4- to OH-.")
    print("The reaction is: Al(OH)3(s) + OH-(aq) <=> Al(OH)4-(aq)")
    print(f"The constant K_c = K_sp * K_f = ({K_sp:.2e}) * ({K_f:.2e}) = {K_c:.2e}")

    # --- Step 3: Set up the charge balance equation ---
    print("\nStep 3: Use the charge balance equation.")
    print("3[Al^3+] + [H+] = [OH-] + [Al(OH)4-]")
    print("Substituting expressions in terms of [OH-]:")
    print(f"3*(K_sp / [OH-]^3) + (K_w / [OH-]) = [OH-] + (K_c * [OH-])")
    
    # --- Step 4: Rearrange into a polynomial equation ---
    print("\nStep 4: Rearrange into a solvable polynomial equation for [OH-].")
    print("Let x = [OH-]. The equation becomes:")
    print(f"(1 + K_c) * x^4 - K_w * x^2 - 3 * K_sp = 0")
    print("This is a quadratic equation in x^2. Let y = x^2.")
    print("A*y^2 + B*y + C = 0, where:")
    
    A = 1 + K_c
    B = -K_w
    C = -3 * K_sp
    print(f"  A = 1 + K_c = {A:.4e}")
    print(f"  B = -K_w = {B:.4e}")
    print(f"  C = -3 * K_sp = {C:.4e}")

    # --- Step 5: Solve the quadratic equation for y = [OH-]^2 ---
    print("\nStep 5: Solve for y = [OH-]^2 using the quadratic formula.")
    discriminant = math.sqrt(B**2 - 4*A*C)
    # We take the positive root for y since y = [OH-]^2 must be positive.
    y = (-B + discriminant) / (2 * A)
    print(f"y = [OH-]^2 = {y:.4e}")
    
    # --- Step 6: Calculate [OH-] ---
    hydroxide_ion_conc = math.sqrt(y)
    print(f"[OH-] = sqrt(y) = {hydroxide_ion_conc:.4e} mol L^-1")

    # --- Step 7: Calculate concentrations of Al species and total solubility ---
    print("\nStep 6: Calculate the concentrations of dissolved aluminum species.")
    al_ion_conc = K_sp / (hydroxide_ion_conc**3)
    complex_ion_conc = K_c * hydroxide_ion_conc
    solubility = al_ion_conc + complex_ion_conc

    print(f"[Al^3+] = K_sp / [OH-]^3 = {al_ion_conc:.4e} mol L^-1")
    print(f"[Al(OH)4-] = K_c * [OH-] = {complex_ion_conc:.4e} mol L^-1")
    
    print("\nStep 7: Calculate the total solubility (s).")
    print("s = [Al^3+] + [Al(OH)4-]")
    print(f"s = {al_ion_conc:.4e} + {complex_ion_conc:.4e}")
    print(f"s = {solubility:.4e} mol L^-1")

    # --- Step 8: Final Answer ---
    print("\nFinal Answer (rounded to two significant figures):")
    print(f"The solubility of Al(OH)3 in pure water is {solubility:.2e} mol L^-1.")
    
    return solubility

solubility_value = calculate_solubility()
# The final answer format as requested
print(f"<<<{solubility_value:.2e}>>>")
