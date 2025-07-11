import sympy

def analyze_soliton_stability():
    """
    Analyzes the stability of a localized soliton in 3D using a scaling argument.
    The Hamiltonian density is H = A*(nabla m)^2 + D*m.(nabla x m).
    """

    # --- Setup ---
    # Define symbols for the scaling analysis.
    # lambda is the dimensionless scaling factor.
    lambda_scale = sympy.symbols('lambda', positive=True)
    # E_ex(1) and E_DMI(1) are the energy values for the unscaled solution (lambda=1).
    # For any non-uniform soliton, the exchange energy E_ex(1) must be positive.
    E_ex_1 = sympy.Symbol('E_ex(1)', positive=True)
    # The sign of the DMI energy E_DMI(1) is not assumed.
    E_dmi_1 = sympy.Symbol('E_DMI(1)')
    # The spatial dimension is 3.
    dim = 3

    print("--- Soliton Stability Analysis using Scaling Argument ---")
    print(f"We are in {dim} spatial dimensions.\n")

    # --- Energy Scaling ---
    # The Heisenberg exchange term ~ (nabla m)^2 scales as lambda^(D-2).
    # The Dzyaloshinskii-Moriya (DMI) term ~ m.(nabla x m) scales as lambda^(D-1).
    E_lambda = E_ex_1 * lambda_scale**(dim - 2) + E_dmi_1 * lambda_scale**(dim - 1)

    print("Step 1: Express total energy E as a function of the scaling factor lambda.")
    print("E(lambda) = E_ex(1) * lambda^(D-2) + E_DMI(1) * lambda^(D-1)")
    print(f"For D={dim}, this becomes:")
    # Use sympy.Eq for pretty printing the equation
    display_eq = sympy.Eq(sympy.Symbol("E(lambda)"), E_lambda)
    print(display_eq)
    print("-" * 50)

    # --- Stationary Condition (First Derivative) ---
    print("Step 2: Find the condition for a stationary solution.")
    print("This requires the first derivative dE/d(lambda) to be zero at lambda=1.\n")
    dE_dlambda = sympy.diff(E_lambda, lambda_scale)
    
    print(f"The first derivative is: dE/d(lambda) = {dE_dlambda}")

    # Evaluate at lambda = 1 to find the condition for an equilibrium solution.
    stationary_condition_eq = dE_dlambda.subs(lambda_scale, 1)
    print("At lambda=1, the condition is:")
    final_stationary_eq_display = sympy.Eq(stationary_condition_eq, 0)
    print(final_stationary_eq_display)

    # To satisfy the user's specific instruction, we output the numbers from this equation.
    # The equation is E_ex(1) + 2*E_DMI(1) = 0
    print("\nThe numerical coefficients in the final stationary equation are:")
    c_ex = stationary_condition_eq.coeff(E_ex_1)
    c_dmi = stationary_condition_eq.coeff(E_dmi_1)
    rhs = 0
    print(f"Coefficient of E_ex(1): {int(c_ex)}")
    print(f"Coefficient of E_DMI(1): {int(c_dmi)}")
    print(f"Right-hand side: {rhs}")
    print("\nThis condition implies E_ex(1) = -2 * E_DMI(1). Since E_ex(1) must be positive,")
    print("it follows that E_DMI(1) must be negative for any stationary soliton.")
    print("-" * 50)

    # --- Stability Condition (Second Derivative) ---
    print("Step 3: Check the condition for stability.")
    print("This requires the second derivative d^2E/d(lambda)^2 to be positive at lambda=1.\n")
    d2E_dlambda2 = sympy.diff(dE_dlambda, lambda_scale)

    print(f"The second derivative is: d^2E/d(lambda)^2 = {d2E_dlambda2}")

    # Evaluate at lambda=1 to find the stability condition.
    stability_value = d2E_dlambda2.subs(lambda_scale, 1)
    print("At lambda=1, the stability condition is:")
    print(f"{stability_value} > 0")
    print(f"This simplifies to: E_DMI(1) > 0.")
    print("-" * 50)
    
    # --- Conclusion ---
    print("Step 4: Combine the conditions and conclude.")
    print("We have found two contradictory requirements for a stable soliton:")
    print("1. From the stationary condition: E_DMI(1) must be NEGATIVE.")
    print("2. From the stability condition: E_DMI(1) must be POSITIVE.")
    
    print("\nThese two conditions cannot be satisfied at the same time.")
    print("\nTherefore, any stationary solution that might exist is not a true energy minimum")
    print("but a saddle point, making it inherently unstable to expansion or collapse.")
    print("\nFinal Answer: It is NOT possible to stabilize localized solitons with this Hamiltonian.")

if __name__ == '__main__':
    analyze_soliton_stability()