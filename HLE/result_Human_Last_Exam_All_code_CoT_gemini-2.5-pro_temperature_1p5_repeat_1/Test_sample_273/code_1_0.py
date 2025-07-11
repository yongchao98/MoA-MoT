import sympy

def analyze_soliton_stability():
    """
    Performs a symbolic scaling analysis (Derrick's theorem) to determine
    the stability of a 3D localized soliton in a system with only
    Heisenberg exchange and DMI interactions.
    """
    # Define symbolic variables for the analysis
    # lmbda: The positive, real scaling parameter for the soliton's size.
    # E_A: The exchange energy component, which must be positive for a non-uniform soliton.
    # E_D: The Dzyaloshinskii-Moriya (DMI) energy component.
    lmbda = sympy.Symbol('lambda', real=True, positive=True)
    E_A = sympy.Symbol('E_A', real=True, positive=True)
    E_D = sympy.Symbol('E_D', real=True)

    # In 3D, the exchange energy scales as lambda^(-1) and the DMI energy as lambda^(-2).
    # The total energy of the scaled configuration is:
    E_scaled = E_A / lmbda + E_D / (lmbda**2)

    # For a stationary soliton solution, the first derivative of the energy with respect to
    # the scaling parameter must be zero at lambda = 1.
    dE_dlmbda = sympy.diff(E_scaled, lmbda)
    stationarity_condition = sympy.Eq(dE_dlmbda.subs(lmbda, 1), 0)

    # Solve the stationarity condition to find the required relation between E_A and E_D.
    # This relationship is a form of the virial theorem for this system.
    virial_relation = sympy.solve(stationarity_condition, E_D)[0]

    # Stability is determined by the sign of the second derivative of the energy.
    # A positive value implies a stable minimum, while a negative value implies an unstable maximum.
    d2E_dlmbda2 = sympy.diff(dE_dlmbda, lmbda)
    stability_at_one = d2E_dlmbda2.subs(lmbda, 1)

    # Substitute the virial relation into the stability expression to get the final result.
    final_stability_result = stability_at_one.subs(E_D, virial_relation)

    # Print the steps of the derivation and the final conclusion.
    print("Derrick's Theorem Scaling Analysis for 3D Soliton Stability")
    print("-" * 60)
    print(f"Energy of scaled configuration: E(lambda) = {E_scaled}")
    print(f"First derivative dE/d(lambda): {dE_dlmbda}")
    print(f"Stationarity condition (dE/d(lambda) = 0 at lambda=1): {stationarity_condition}")
    print(f"This implies a virial relation between energy terms: E_D = {virial_relation}")
    print("-" * 60)
    print(f"Second derivative d^2E/d(lambda)^2: {d2E_dlmbda2}")
    print(f"Value at lambda=1 (before substitution): {stability_at_one}")
    print("Substituting the virial relation gives the stability criterion.")

    final_eq_lhs = "d^2E/d(lambda)^2 |_(lambda=1)"
    final_eq_rhs = final_stability_result

    # Fulfilling the request to "print each number in the final equation"
    # by breaking the symbolic expression into its coefficient and variable parts.
    coeff, var = final_eq_rhs.as_coeff_Mul()

    print(f"\nThe final stability equation is: {final_eq_lhs} = {final_eq_rhs}")
    print("\nBreaking down the terms of the final equation's result:")
    print(f"Coefficient: {coeff}")
    print(f"Variable Part: {var}")

    print("\nConclusion:")
    print("Since E_A (the exchange energy) must be positive for any non-uniform soliton,")
    print("the second derivative (-1 * E_A) is always negative.")
    print("This indicates an instability (an energy maximum), proving that a localized")
    print("soliton cannot be stabilized in this 3D model.")

if __name__ == '__main__':
    analyze_soliton_stability()