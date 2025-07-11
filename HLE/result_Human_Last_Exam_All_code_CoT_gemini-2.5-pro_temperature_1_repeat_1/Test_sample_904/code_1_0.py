import sympy

def solve_fluid_equation():
    """
    Derives the governing linear equation for the interfacial shape xi(r)
    and identifies the coefficients A(r) and B(r).
    """
    # Define the symbolic variables and functions
    r = sympy.Symbol('r', positive=True, real=True)
    gamma = sympy.Symbol('gamma', positive=True, real=True)
    xi = sympy.Function('xi')(r)
    P_el = sympy.Function('P_el')(r)

    print("Step 1: State the Young-Laplace equation: Delta_P = gamma * kappa")
    print(f"Here, Delta_P is the electrostatic pressure P_el(r) and gamma is surface tension.")

    # Define the linearized mean curvature kappa for an axisymmetric surface xi(r)
    kappa = xi.diff(r, 2) + (1/r) * xi.diff(r, 1)
    print("\nStep 2: Define the linearized mean curvature kappa in cylindrical coordinates:")
    print(f"kappa = {kappa}")

    # Formulate the governing differential equation based on Young-Laplace
    # gamma * kappa - P_el = 0
    governing_eq = sympy.Eq(gamma * kappa - P_el, 0)
    print("\nStep 3: Formulate the governing equation:")
    print(governing_eq)

    # Rearrange the equation into the form A(r)*xi'' + B(r)*xi' + C = 0
    # To achieve a standard form, we multiply by r to clear the denominator.
    canonical_eq = sympy.simplify(r * governing_eq.lhs)
    print("\nStep 4: Multiply by r to clear the denominator:")
    print(sympy.Eq(canonical_eq, 0))

    # Divide by the constant gamma to get the most standard form.
    # The influence of gamma is now in the C(r) term.
    final_eq = sympy.expand(canonical_eq / gamma)
    print("\nStep 5: Divide by the constant gamma to obtain the canonical form:")
    print(sympy.Eq(final_eq, 0))
    print("\nThis equation is in the form A(r)*xi'' + B(r)*xi' + C(r) = 0.")

    # Extract the coefficients A(r) and B(r) from the final equation
    xi_dd = xi.diff(r, 2)
    xi_d = xi.diff(r, 1)

    A_r = final_eq.coeff(xi_dd)
    B_r = final_eq.coeff(xi_d)

    print("\nStep 6: Identify the coefficients A(r) and B(r).")
    print("\n--- Final Result ---")
    print(f"The coefficient A(r) is the term multiplying the second derivative d^2xi/dr^2.")
    print(f"A(r) = {A_r}")
    print("\nThe coefficient B(r) is the term multiplying the first derivative dxi/dr.")
    print(f"B(r) = {B_r}")

solve_fluid_equation()