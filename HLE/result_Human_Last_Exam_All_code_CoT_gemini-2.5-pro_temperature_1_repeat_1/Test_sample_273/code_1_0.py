import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a 3D soliton using a scaling argument (Derrick's Theorem).

    The energy functional is E = integral[ A*(grad m)^2 + D*m.curl(m) ] dV.
    We test for stability by considering a scaled solution m_lambda(x) = m_0(lambda*x),
    which effectively changes the size of the soliton. The energy of this scaled
    solution, E(lambda), can be shown to have the following form based on how the
    gradient and volume elements transform under this scaling in 3 dimensions:
    E(lambda) = E_A * lambda^(-1) + E_D * lambda^(-2)
    where E_A is the total exchange energy and E_D is the total DMI energy of the
    original (lambda=1) solution.
    """
    # Define symbolic variables
    lmbda = sympy.Symbol('lambda', positive=True)
    E_A = sympy.Symbol('E_A', positive=True) # Exchange energy, must be positive for a non-uniform texture
    E_D = sympy.Symbol('E_D', real=True)      # DMI energy

    # Define the energy of the scaled solution
    E = E_A * lmbda**(-1) + E_D * lmbda**(-2)

    print("This script analyzes the stability of a 3D magnetic soliton based on the provided Hamiltonian.")
    print("We use a scaling argument (Derrick's Theorem).\n")
    print(f"The energy of a scaled configuration is given by: E(lambda) = {E}\n")

    # First derivative for finding stationary points (extrema)
    dE_dlmbda = sympy.diff(E, lmbda)
    print(f"The first derivative of energy with respect to the scaling parameter lambda is:")
    print(f"dE/dlambda = {dE_dlmbda}\n")

    # For a solution to be a stationary point, dE/dlambda must be zero at lambda=1
    stationary_condition_eq = sympy.Eq(dE_dlmbda.subs(lmbda, 1), 0)
    print("For a stationary solution (a candidate for a stable soliton), dE/dlambda must be 0 at lambda = 1.")
    print(f"This gives the condition (virial theorem): {stationary_condition_eq}\n")

    # Solve for the relationship between E_A and E_D
    virial_relation = sympy.solve(stationary_condition_eq, E_A)[0]
    print(f"From this condition, we find the relation: E_A = {virial_relation}")
    print("Since E_A (from the exchange term) must be positive for any non-uniform texture, E_D must be negative.\n")

    # Second derivative for testing stability (minimum vs maximum)
    d2E_dlmbda2 = sympy.diff(dE_dlmbda, lmbda)
    print("To test for stability, we examine the sign of the second derivative at lambda = 1.")
    print(f"The second derivative is: d^2E/dlambda^2 = {d2E_dlmbda2}\n")

    # Evaluate the second derivative at lambda=1
    second_derivative_at_1 = d2E_dlmbda2.subs(lmbda, 1)
    print(f"At lambda = 1, the second derivative evaluates to: {second_derivative_at_1}\n")

    # Substitute the virial theorem into the second derivative to determine its sign
    print("Now, we substitute the stationary point condition (E_A = -2*E_D) into the second derivative:")
    # Manually print the steps for clarity
    print(f"Expression: 2*E_A + 6*E_D")
    print(f"Substitute E_A: 2*({virial_relation}) + 6*E_D")
    final_expr = 2*virial_relation + 6*E_D
    print(f"Final Result: {final_expr}\n")

    print("--- Conclusion ---")
    print("For a stable soliton, the energy must be at a local minimum, meaning the second derivative must be positive.")
    print(f"This would require the final expression, {final_expr}, to be positive, which implies E_D > 0.")
    print("However, the condition for the existence of a stationary point (E_A = -2*E_D) requires E_D to be negative, because E_A must be positive.")
    print("These two conditions are contradictory.")
    print("The stationary point is an energy maximum with respect to scaling (since 2*E_D < 0), not an energy minimum.")
    print("Therefore, a localized soliton described by this Hamiltonian is fundamentally unstable.")

# Run the analysis
analyze_soliton_stability()