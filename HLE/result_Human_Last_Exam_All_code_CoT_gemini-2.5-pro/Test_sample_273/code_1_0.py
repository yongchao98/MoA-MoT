import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in 3D using a scaling argument.

    The Hamiltonian density is H = A(∇m)² + Dm·(∇×m).
    We analyze the stability of a supposed localized solution m(r) by considering
    a scaled version m(r/λ) and examining the total energy E(λ).
    """

    # Define symbolic variables
    # lmbda is the scaling parameter
    # E_A and E_D are the energy values of the unscaled solution (λ=1)
    # for the exchange and DMI terms, respectively.
    lmbda = sympy.symbols('lambda', real=True, positive=True)
    E_A = sympy.Symbol('E_A', real=True)
    E_D = sympy.Symbol('E_D', real=True)

    print("Step 1: Define the total energy of the scaled solution E(λ).")
    print("In 3D, the exchange energy E_A scales as λ¹.")
    print("In 3D, the DMI energy E_D scales as λ².")
    
    # Total energy of the scaled solution
    E_lambda = E_A * lmbda + E_D * lmbda**2
    print(f"\nThe total energy is E(λ) = {E_lambda}\n")

    print("Step 2: Find the equilibrium condition.")
    print("For a stationary solution, the first derivative dE/dλ must be zero at λ=1.")
    
    # First derivative of E with respect to lambda
    dE_dlmbda = sympy.diff(E_lambda, lmbda)
    print(f"The first derivative is dE/dλ = {dE_dlmbda}")

    # Set the first derivative to zero at lambda=1
    equilibrium_condition = sympy.Eq(dE_dlmbda.subs(lmbda, 1), 0)
    print(f"At λ=1, the equilibrium condition is: {equilibrium_condition}")
    print("This is the virial theorem for this system.\n")

    print("Step 3: Find the stability condition.")
    print("For the solution to be a stable minimum, the second derivative d²E/dλ² must be positive at λ=1.")

    # Second derivative of E with respect to lambda
    d2E_dlmbda2 = sympy.diff(dE_dlmbda, lmbda)
    print(f"The second derivative is d²E/dλ² = {d2E_dlmbda2}")

    # Evaluate at lambda=1
    stability_value = d2E_dlmbda2.subs(lmbda, 1)
    print(f"At λ=1, the second derivative is {stability_value}.")
    print(f"Thus, for stability, we must have {stability_value} > 0, which implies E_D > 0.\n")

    print("Step 4: Check for contradiction.")
    print("Let's analyze the two conditions together.")
    print(f"Condition 1 (Equilibrium): {equilibrium_condition}")
    print("Condition 2 (Stability): E_D > 0\n")

    print("From Condition 1, we can express E_A in terms of E_D:")
    # Solve for E_A from the equilibrium condition
    E_A_relation = sympy.solve(equilibrium_condition, E_A)[0]
    print(f"E_A = {E_A_relation}")

    print("\nNow, we must consider the physical nature of the exchange energy E_A.")
    print("The exchange term is ∫A(∇m)²dV, where A > 0. For any non-uniform soliton, (∇m)² > 0 somewhere.")
    print("Therefore, the total exchange energy E_A must be positive: E_A > 0.")

    print("\nIf E_A > 0, then from the relation E_A = -2*E_D, we must have:")
    print("-2*E_D > 0, which implies E_D < 0.")

    print("\nConclusion:")
    print("The stability condition requires E_D > 0.")
    print("However, the equilibrium condition combined with the fact that E_A > 0 forces E_D < 0.")
    print("These two conditions are contradictory.")
    print("\nTherefore, no stable localized soliton can exist in this 3D model.")

if __name__ == '__main__':
    analyze_soliton_stability()