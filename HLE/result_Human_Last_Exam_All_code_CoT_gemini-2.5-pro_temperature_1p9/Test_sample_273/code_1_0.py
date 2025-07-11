import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in 3D using the Hobart-Derrick theorem
    for a Hamiltonian with only exchange and DMI terms.
    """

    # --- 1. Define symbolic variables ---
    # λ is the scaling factor
    # E_ex and E_dmi are the energy values of the unscaled solution (λ=1)
    # D is the number of spatial dimensions
    L = sympy.Symbol('λ', positive=True)
    E_ex = sympy.Symbol('E_ex', positive=True) # Exchange energy must be positive for non-trivial texture
    E_dmi = sympy.Symbol('E_dmi', real=True)   # DMI energy can be positive or negative
    D = 3

    # --- 2. State the problem and scaling rules ---
    print("Analyzing the stability of a 3D localized soliton.")
    print("The Hamiltonian is H = integral[A*(grad m)^2 + D*m.curl(m)] dV")
    print("-" * 60)
    print("According to the Hobart-Derrick theorem, we analyze the energy under a scaling transformation.")
    print(f"In {D} dimensions:")
    exponent_ex = D - 2
    exponent_dmi = D - 1
    print(f"  - The exchange term E_ex scales as λ^({D}-2) = λ^{exponent_ex}.")
    print(f"  - The DMI term E_dmi scales as λ^({D}-1) = λ^{exponent_dmi}.")
    print("-" * 60)

    # --- 3. Formulate the scaled total energy E(λ) ---
    E_L = E_ex * L**exponent_ex + E_dmi * L**exponent_dmi

    print("The scaled total energy E(λ) is:")
    print(f"E(λ) = E_ex * λ^{exponent_ex} + E_dmi * λ^{exponent_dmi}")
    print("-" * 60)

    # --- 4. First condition: stationary solution (dE/dλ at λ=1 must be 0) ---
    dE_dL = sympy.diff(E_L, L)
    print("First derivative of energy, dE/dλ:")
    print(f"dE/dλ = {dE_dL}")

    # Evaluate at L=1
    dE_dL_at_1 = dE_dL.subs(L, 1)
    print("\nFor a stationary solution, dE/dλ at λ=1 must be zero:")
    print(f"Condition 1: {dE_dL_at_1} = 0")
    
    # This gives a relationship between E_ex and E_dmi
    print(f"\nThis gives the equilibrium equation:")
    # Printing each number in the final equation as requested
    print(f"({exponent_ex}) * E_ex + ({exponent_dmi}) * E_dmi = 0")

    # --- 5. Second condition: stability (d^2E/dλ^2 at λ=1 must be >= 0) ---
    d2E_dL2 = sympy.diff(dE_dL, L)
    print("-" * 60)
    print("Second derivative of energy, d^2E/dλ^2:")
    print(f"d^2E/dλ^2 = {d2E_dL2}")

    # Evaluate at L=1
    d2E_dL2_at_1 = d2E_dL2.subs(L, 1)
    print("\nFor a stable solution, d^2E/dλ^2 at λ=1 must be non-negative:")
    print(f"Condition 2: {d2E_dL2_at_1} >= 0")
    print(f"This implies that ({exponent_dmi} * ({exponent_dmi} - 1)) * E_dmi >= 0, which simplifies to {2*E_dmi} >= 0.")
    print("-" * 60)
    
    # --- 6. Combine conditions and conclude ---
    print("Let's analyze the conditions together:")
    print(f"Condition 1:  E_ex + 2 * E_dmi = 0  =>  E_ex = -2 * E_dmi")
    print(f"Condition 2:  2 * E_dmi >= 0         =>  E_dmi >= 0")

    print("\nPhysical constraints:")
    print("For any localized soliton (a non-uniform texture), the exchange energy E_ex = integral[A*(grad m)^2] dV must be strictly positive (E_ex > 0).")
    
    print("\nChecking for contradictions:")
    print("If E_dmi >= 0 (from Condition 2), then for Condition 1 to hold, E_ex = -2 * E_dmi must be non-positive (E_ex <= 0).")
    print("This contradicts the physical requirement that E_ex must be strictly positive for a soliton.")
    print("The only way to satisfy both requirements is if E_ex=0 and E_dmi=0, which corresponds to the uniform ground state, not a localized soliton.")
    
    print("\nFinal conclusion:")
    print("No. A localized soliton cannot be stabilized in 3D with only Heisenberg exchange and DMI terms, as it is unstable against collapse.")

if __name__ == '__main__':
    analyze_soliton_stability()