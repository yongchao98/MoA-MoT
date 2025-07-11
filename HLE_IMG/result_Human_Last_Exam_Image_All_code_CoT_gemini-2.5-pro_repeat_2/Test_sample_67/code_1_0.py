import sympy

def solve_electron_energy_problem():
    """
    This function derives the minimum energy for an electron to cause impact ionization
    based on the provided energy band structure.
    """

    # --- Symbolic Setup ---
    # Using strings for symbolic representation in the explanation
    E_g_str = "E_g"   # Band gap energy
    m_star_str = "m*" # Effective mass
    h_bar_str = "ħ"   # Reduced Planck constant
    k_i, k_1, k_2, k_h = sympy.symbols('k_i k_1 k_2 k_h') # Wave vectors

    print("Step 1: Identify the physical process and conservation laws.")
    print("The process is impact ionization. An initial electron (e_i) promotes an electron from the valence band,")
    print("resulting in a final state of two electrons (e_1, e_2) in the conduction band and one hole (h) in the valence band.")
    print("Process: e_i -> e_1 + e_2 + h\n")
    print("We apply conservation of energy and momentum (wave vector k).")
    print("Let K(k) = (ħ²k²)/(2m*) be the kinetic energy term.\n")

    print("Step 2: Formulate the conservation equations.")
    print("Energy Conservation: E_i = E_1 + E_2 + E_h")
    print(f"E_i is the initial energy of electron 1. E_1 and E_2 are final electron energies. E_h is the hole energy.")
    print(f"( {E_g_str} + K(k_i) ) = ( {E_g_str} + K(k_1) ) + ( {E_g_str} + K(k_2) ) + K(k_h)")
    print(f"This simplifies to: K(k_i) = {E_g_str} + K(k_1) + K(k_2) + K(k_h)  (Equation 1)\n")

    print("Momentum Conservation:")
    print("k_i = k_1 + k_2 - k_h  (Equation 2)")
    print("(The momentum of the created hole k_h is subtracted from the final electron momenta)\n")

    print("Step 3: Solve for the minimum initial energy E_i.")
    print("The initial energy is E_i = E_g + K(k_i). To minimize E_i, we must minimize K(k_i).")
    print("From Eq. 1, this means we must minimize the sum of the kinetic energies of the final state particles: K_sum = K(k_1) + K(k_2) + K(k_h).")
    print("This minimization is subject to the constraint imposed by both conservation laws.\n")
    
    print("Step 4: Minimize final state energy subject to the constraint.")
    print("Combining the conservation laws leads to a constraint on the final state wave vectors:")
    print(f"k_1·k_2 - k_1·k_h - k_2·k_h = ({m_star_str}*{E_g_str}) / {h_bar_str}²")
    print("To minimize the sum of kinetic energies K_sum under this constraint, all final particles must be collinear.")
    print("The optimal configuration that requires the least energy is found to be k_1 = k, k_2 = k, and k_h = -k.\n")

    print("Step 5: Calculate the final result.")
    print("Substituting this configuration into the constraint gives: k² - k(-k) - k(-k) = 3k².")
    print(f"So, 3k² = ({m_star_str}*{E_g_str}) / {h_bar_str}².\n")

    print("The minimum sum of final kinetic energies is:")
    print(f"K_sum_min = K(k) + K(k) + K(-k) = ({h_bar_str}²/(2*{m_star_str})) * (k² + k² + (-k)²) = ({h_bar_str}²/(2*{m_star_str})) * (3k²)")
    print(f"Substituting the result for 3k²:")
    print(f"K_sum_min = ({h_bar_str}²/(2*{m_star_str})) * (({m_star_str}*{E_g_str})/{h_bar_str}²) = (1/2)*{E_g_str}\n")
    
    # Perform the final calculation
    E_kin_final_coeff = 0.5
    E_i_min_coeff = 1 + E_kin_final_coeff

    print("Finally, the minimum required initial energy for electron 1 is found from Equation 1:")
    print(f"E_i_min = {E_g_str} + K_sum_min")
    print(f"E_min = {1} * {E_g_str} + {E_kin_final_coeff} * {E_g_str} = {E_i_min_coeff} * {E_g_str}")

if __name__ == '__main__':
    solve_electron_energy_problem()