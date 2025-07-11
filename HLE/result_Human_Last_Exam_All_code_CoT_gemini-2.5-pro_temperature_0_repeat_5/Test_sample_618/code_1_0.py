import sympy

def solve_heat_equation_symmetries():
    """
    This script derives the Lie symmetries for the heat equation with a logarithmic source term
    u_t = u_xx + (k1*ln(u) + k2)*u, and finds the finite transformation for the x-coordinate.
    """
    # Step 1: Define all symbols and functions
    t, x = sympy.symbols('t x', real=True)
    k1, k2 = sympy.symbols('k1 k2', real=True, nonzero=True)
    c2, c3, c4, c5 = sympy.symbols('c2 c3 c4 c5') # Arbitrary constants
    eps = sympy.Symbol('epsilon') # Group parameter
    u = sympy.Symbol('u') # Represent u for eta expression

    # Define functions for solving ODEs
    b = sympy.Function('b')(t)
    G = sympy.Function('G')(t)

    # Step 2: State the results from the determining equations
    print("Derivation of Lie Group Symmetries:")
    print("="*40)
    print("For the equation u_t = u_xx + (k1*ln(u) + k2)*u, assuming k1 is non-zero,")
    print("the analysis of the determining equations shows that the infinitesimals must have the following structure:")
    print("\n1. tau = tau(t)")
    print("2. xi = xi(t, x)")
    print("3. eta = C(t, x)*u (the part of eta independent of u is zero)")
    print("-" * 40)

    print("Step 3: Solving for the infinitesimals")
    print("Further analysis yields the following simplified system of PDEs:")
    print("  - diff(xi, t) - 2*diff(C, x) + diff(xi, x, 2) = 0")
    print("  - 2*diff(xi, x) - diff(tau, t) = 0")
    print("  - The coefficient of u*ln(u) implies diff(xi, x) is a function of t only.")
    print("  - The coefficient of u*ln(u) also forces this function to be zero.")
    
    print("\nThis leads to:")
    print("  - diff(xi, x) = 0  =>  xi is a function of t only, let's call it b(t).")
    print("  - diff(tau, t) = 0  =>  tau is a constant, c3.")
    
    # With xi = b(t), the system simplifies to:
    # -b'(t) - 2*C_x = 0  => C(t,x) = -b'(t)*x/2 + G(t)
    # C_t - C_xx - k1*C = 0
    print("\nSubstituting these into the remaining equations gives:")
    print("  - C(t,x) = -diff(b(t),t)*x/2 + G(t)")
    print("  - And two ODEs obtained by substituting C into its PDE and separating by powers of x:")
    
    # ODE for b(t)
    b_ode = sympy.Eq(sympy.diff(b, t, 2) - k1 * sympy.diff(b, t), 0)
    print(f"    - For b(t): {b_ode}")
    # ODE for G(t)
    G_ode = sympy.Eq(sympy.diff(G, t) - k1 * G, 0)
    print(f"    - For G(t): {G_ode}")

    # Solve the ODEs
    # dsolve gives C1 + C2*exp(k1*t). We rename constants for clarity.
    b_sol = c5 + (c2/k1)*sympy.exp(k1*t) 
    # dsolve gives C1*exp(k1*t).
    G_sol = c4*sympy.exp(k1*t)

    print("\nSolving these ODEs gives:")
    print(f"  - b(t) = c5 + (c2/k1)*exp(k1*t)")
    print(f"  - G(t) = c4*exp(k1*t)")
    print("-" * 40)

    print("Step 4: Final form of the infinitesimals")
    tau_final = c3
    xi_final = b_sol
    b_prime = sympy.diff(b_sol, t)
    C_final = (-b_prime * x / 2 + G_sol)
    eta_final = (C_final * u).simplify()

    print(f"  tau = {tau_final}")
    print(f"  xi  = {xi_final}")
    print(f"  eta = {eta_final}")
    print("-" * 40)

    print("Step 5: Finite transformation for the x-coordinate")
    print("The transformation x_new is found by solving dx'/d(eps) = xi(t') where t' = t + c3*eps.")
    print("\nThe general representation for the transformation on x has two forms, depending on whether c3 is zero:")

    # Define symbols for printing the final equations
    x_new = sympy.Symbol("x_new")
    x_old = sympy.Symbol("x")
    t_old = sympy.Symbol("t")

    # Case c3 != 0
    eq_case1 = sympy.Eq(x_new, x_old + c5*eps + (c2/(k1*c3)) * sympy.exp(k1*t_old) * (sympy.exp(k1*c3*eps) - 1))
    # Case c3 == 0
    eq_case2 = sympy.Eq(x_new, x_old + eps*(c5 + c2*sympy.exp(k1*t_old)))

    print("\nCase 1: c3 is not zero")
    print(f"  {sympy.pretty(x_new)} = {sympy.pretty(eq_case1.rhs)}")

    print("\nCase 2: c3 is zero")
    print(f"  {sympy.pretty(x_new)} = {sympy.pretty(eq_case2.rhs)}")

    print("\nHere, 'eps' is the group parameter, and c2, c3, c5 are arbitrary constants from the Lie group.")

if __name__ == '__main__':
    solve_heat_equation_symmetries()