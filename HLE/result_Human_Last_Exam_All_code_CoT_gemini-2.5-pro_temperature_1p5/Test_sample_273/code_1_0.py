import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in 3D using a scaling argument
    (a form of the Hobart-Derrick theorem).

    The total energy of a magnetic texture is given by the Hamiltonian:
    H = ∫[A(∇m)² + Dm·(∇×m)] dV

    If we consider a localized texture and scale its size by a factor λ, the
    total energy E as a function of λ can be written as:
    E(λ) = E_A * λ + E_D * λ²

    where:
    - E_A is the exchange energy component (> 0 for a non-uniform texture).
    - E_D is the DMI component.
    """

    # Use the sympy library for symbolic mathematics to make the logic transparent.
    l, E_A_sym, E_D_sym = sympy.symbols('lambda E_A E_D')

    # Define the symbolic energy function based on the scaling properties.
    E_sym = E_A_sym * l + E_D_sym * l**2

    # Calculate the first and second derivatives with respect to the size λ.
    dE_dl_sym = sympy.diff(E_sym, l)
    d2E_dl2_sym = sympy.diff(dE_dl_sym, l)

    print("--- Stability Analysis of a 3D Soliton ---")
    print("\nThe energy of a scaled soliton is described by the equation:")
    # The final equation string is built here.
    equation_str = f"E(lambda) = (E_A * lambda) + (E_D * lambda^2)"
    print(equation_str)

    print("\n--- Step 1: Condition for Equilibrium ---")
    print("An equilibrium state (a potential soliton solution) exists if the energy has a stationary point for lambda > 0.")
    print("This occurs when the 'force' (the first derivative of energy) is zero.")
    print(f"dE/dlambda = {dE_dl_sym}")
    print("For dE/dlambda to be zero at a positive lambda, E_A and E_D must have opposite signs.")
    print("Since exchange energy E_A must be positive, the DMI energy E_D must be NEGATIVE (E_D < 0).")


    print("\n--- Step 2: Condition for Stability ---")
    print("For the equilibrium to be stable, it must be an energy MINIMUM.")
    print("This requires the second derivative of energy to be POSITIVE.")
    print(f"d^2E/dlambda^2 = {d2E_dl2_sym}")
    print("For d^2E/dlambda^2 to be positive, E_D must be POSITIVE (E_D > 0).")

    print("\n--- Conclusion from Analysis ---")
    print("We have a contradiction:")
    print("  - Equilibrium requires: E_D < 0")
    print("  - Stability requires:   E_D > 0")
    print("It is impossible to satisfy both conditions at the same time.")
    print("Therefore, no stable localized soliton can exist with only these two energy terms.")

    print("\n--- Numerical Example ---")
    # To demonstrate, we choose numbers that satisfy the equilibrium condition at lambda=1.
    # This means E_A + 2*E_D = 0, so E_A = -2*E_D.
    # We must pick E_D < 0, so let's choose E_D = -1.0.
    E_A_val = 2.0
    E_D_val = -1.0
    print(f"Let's assume a hypothetical equilibrium exists and pick values: E_A = {E_A_val}, E_D = {E_D_val}")

    # Calculate the numerical value of the second derivative.
    d2E_dl2_val = 2 * E_D_val

    print("\nThe final equation for the stability check is the second derivative:")
    print(f"d^2E/dlambda^2 = 2 * E_D")
    print("\nPlugging in the numbers:")
    print(f"d^2E/dlambda^2 = 2 * ({E_D_val}) = {d2E_dl2_val}")
    print(f"\nSince the result is {d2E_dl2_val} (negative), the equilibrium is UNSTABLE.")
    print("\nThe soliton would either collapse to a point or expand indefinitely.")

if __name__ == '__main__':
    analyze_soliton_stability()