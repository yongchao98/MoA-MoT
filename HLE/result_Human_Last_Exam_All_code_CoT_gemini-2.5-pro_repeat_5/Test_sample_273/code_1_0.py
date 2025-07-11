import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a 3D localized soliton using a scaling argument (Derrick's theorem).

    The Hamiltonian density is H = A*(grad(m))^2 + D*m.dot(curl(m)).
    We analyze the stability of a solution by scaling its size by a factor lambda.
    """
    # Define symbolic variables
    # L represents the scaling factor lambda
    # E_ex1 and E_DMI1 are the exchange and DMI energies of the unscaled solution (lambda=1)
    L = sympy.Symbol('L', positive=True)
    E_ex1 = sympy.Symbol('E_ex1', positive=True) # Exchange energy must be positive
    E_DMI1 = sympy.Symbol('E_DMI1')

    # In 3D, the energy terms scale as follows:
    # Exchange energy (two derivatives): scales as L^(3-2) = L^1
    # DMI energy (one derivative): scales as L^(3-1) = L^2
    E_L = L * E_ex1 + L**2 * E_DMI1

    print("--- Soliton Stability Analysis (Derrick's Theorem) ---")
    print(f"Total energy E as a function of scaling factor L (lambda):")
    print(f"E(L) = {E_L}")
    print("-" * 55)

    # First condition for a stable solution: It must be an energy extremum.
    # The first derivative of energy with respect to L must be zero at L=1.
    dE_dL = sympy.diff(E_L, L)
    print("First derivative of energy with respect to L:")
    print(f"dE/dL = {dE_dL}")

    # Evaluate at L=1 to find the stationary condition (virial theorem)
    virial_condition = dE_dL.subs(L, 1)
    print("\nCondition for a stationary point (dE/dL at L=1 = 0):")
    # Extract and print the coefficients from the equation
    c_ex = virial_condition.coeff(E_ex1)
    c_dmi = virial_condition.coeff(E_DMI1)
    print(f"({c_ex})*E_ex1 + ({c_dmi})*E_DMI1 = 0")

    # This condition relates E_ex1 and E_DMI1
    relation = sympy.solve(virial_condition, E_DMI1)
    print(f"This implies that for a stationary solution, E_DMI1 = {relation[0]}")
    print("-" * 55)

    # Second condition for a stable solution: It must be an energy MINIMUM.
    # The second derivative of energy with respect to L must be positive at L=1.
    d2E_dL2 = sympy.diff(dE_dL, L)
    print("Second derivative of energy with respect to L:")
    print(f"d^2E/dL^2 = {d2E_dL2}")

    stability_value = d2E_dL2.subs(L, 1)
    print("\nCondition for stability (d^2E/dL^2 at L=1 > 0):")
    # Extract and print the coefficient from the equation
    c_dmi_2 = stability_value.coeff(E_DMI1)
    print(f"({c_dmi_2})*E_DMI1 > 0  =>  E_DMI1 > 0")
    print("-" * 55)

    # Final conclusion: Check for contradiction
    print("--- Conclusion ---")
    print("1. For any non-uniform soliton, the exchange energy E_ex1 must be positive.")
    print(f"2. The stationary condition (E_ex1 + {c_dmi}*E_DMI1 = 0) requires E_DMI1 = {relation[0]}.")
    print("   Since E_ex1 > 0, this means E_DMI1 must be negative.")
    print("3. The stability condition requires E_DMI1 to be positive.")
    print("\nThe condition for a stationary point (E_DMI1 < 0) contradicts the condition for stability (E_DMI1 > 0).")
    print("\nTherefore, no stable localized soliton can exist with only these two energy terms in 3D.")

if __name__ == '__main__':
    analyze_soliton_stability()