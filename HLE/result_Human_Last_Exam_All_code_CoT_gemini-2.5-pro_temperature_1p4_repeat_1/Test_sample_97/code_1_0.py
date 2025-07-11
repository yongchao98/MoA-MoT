import math

def solve_solubility():
    """
    Calculates the molar solubility of Al(OH)3 in pure water,
    considering both dissolution and complex ion formation.
    """
    # Given constants
    K_sp = 5.3e-27
    K_f = 1.1e31
    K_w = 1.0e-14

    # --- Step 1: Explain the chemistry and define derived constants ---
    print("The solubility of Al(OH)3 in pure water is determined by two main equilibria:")
    print("1. Dissolution: Al(OH)3(s) <=> Al^3+(aq) + 3 OH^-(aq)")
    print(f"   K_sp = [Al^3+][OH^-]^3 = {K_sp:.2e}\n")
    
    # The overall reaction for forming the complex from the solid is useful.
    # Al(OH)3(s) <=> Al^3+(aq) + 3OH^-(aq)      (K_sp)
    # Al^3+(aq) + 4OH^-(aq) <=> Al(OH)4^-(aq)   (K_f)
    # Net: Al(OH)3(s) + OH^-(aq) <=> Al(OH)4^-(aq) (K_c)
    K_c = K_sp * K_f
    print("A second important equilibrium is the dissolution to form the tetrahydroxoaluminate(III) complex:")
    print("Al(OH)3(s) + OH^-(aq) <=> Al(OH)4^-(aq)")
    print(f"   K_c = K_sp * K_f = ({K_sp:.2e}) * ({K_f:.2e}) = {K_c:.2e}\n")

    # --- Step 2: Set up and solve the charge balance equation for [OH-] ---
    print("The total molar solubility (S) is the sum of the dissolved aluminum species:")
    print("S = [Al^3+] + [Al(OH)4^-]\n")
    print("To find the concentrations, we solve for [OH-] using the charge balance equation:")
    print("[H+] + 3[Al^3+] = [OH^-] + [Al(OH)4^-]\n")
    print("Substituting the equilibrium expressions (using K_sp, K_w, K_c) gives a polynomial in [OH-].")
    print("Let y = [OH-]^2, the equation is:")
    print(f"(1 + K_c) * y^2 - K_w * y - 3 * K_sp = 0\n")

    # Coefficients for the quadratic equation ay^2 + by + c = 0
    a = 1 + K_c
    b = -K_w
    c = -3 * K_sp

    # Solve the quadratic equation for y = [OH-]^2 using the quadratic formula
    # We take the positive root because y must be positive.
    discriminant = b**2 - 4*a*c
    y = (-b + math.sqrt(discriminant)) / (2*a)
    
    # Calculate [OH-]
    hydroxide_conc = math.sqrt(y)
    print(f"Solving this equation gives [OH^-] = {hydroxide_conc:.3e} mol L^-1\n")

    # --- Step 3: Calculate concentrations of dissolved species and total solubility ---
    print("Using this hydroxide concentration, we calculate the concentration of each aluminum species:")
    al_ion_conc = K_sp / (hydroxide_conc**3)
    complex_ion_conc = K_c * hydroxide_conc
    total_solubility = al_ion_conc + complex_ion_conc

    # --- Step 4: Print the final equation and answer ---
    print("The final equation for the total solubility (S) is the sum of these concentrations:")
    print(f"S = [Al^3+] + [Al(OH)4^-]")
    print(f"S = {al_ion_conc:.3e} mol L^-1 + {complex_ion_conc:.3e} mol L^-1")
    
    # Final answer formatted to two significant figures, as per the input constants.
    print(f"S = {total_solubility:.2e} mol L^-1")
    
    # Print the answer in the required format
    print(f"<<<{total_solubility:.2g}>>>")

if __name__ == "__main__":
    solve_solubility()