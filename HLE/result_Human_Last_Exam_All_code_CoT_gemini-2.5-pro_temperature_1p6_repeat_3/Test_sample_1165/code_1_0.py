import sympy as sp

def solve_ode_fluctuations():
    """
    This script calculates the maximum magnitude of fluctuations R as a function of epsilon.
    It follows the plan outlined above, using symbolic mathematics.
    """
    # Step 1: Define symbolic variables
    x, xi, L, epsilon = sp.symbols('x xi L epsilon', real=True, positive=True)
    
    print("Step 1: Define the problem setup.")
    print(f"The domain is [0, L], where L = 1/epsilon.")
    print(f"The leading-order solution is y0(x) = 1 - epsilon*x.")
    print(f"The fluctuation v(x) solves v'' - epsilon*v' = epsilon^2 * (Sum(delta(x-zi)) - rho).")
    print("We approximate this with v'' = epsilon^2 * (Sum(delta(x-zi)) - rho).\n")

    # Step 2: Define the Green's function G0 for the operator d^2/dx^2
    # with boundary conditions G0(0,xi)=0 and G0(L,xi)=0.
    G0 = sp.Piecewise(
        (xi * (x - L) / L, x > xi),
        (x * (xi - L) / L, x < xi)
    )
    
    print("Step 2: Define the Green's function G0(x, xi) for the d^2/dx^2 operator.")
    print(f"G0(x, xi) = {G0}\n")
    
    # Step 3: Calculate the variance of G0(x, z) where z is a random variable
    # uniformly distributed on [0, L]. This is needed for Var[v(x)].
    # We assume z_i are i.i.d. with density p(xi) = 1/L.
    # Var[G0(x,z)] = E[G0(x,z)^2] - (E[G0(x,z)])^2
    
    E_G0_sq = sp.integrate(G0**2, (xi, 0, L)) / L
    E_G0 = sp.integrate(G0, (xi, 0, L)) / L
    
    Var_G0 = sp.simplify(E_G0_sq - E_G0**2)
    
    print("Step 3: Calculate the variance of G0(x, z) for a single random point z.")
    print(f"Var[G0(x, z)] = E[G0^2] - (E[G0])^2 = {Var_G0}\n")
    
    # Step 4: Calculate the variance of the fluctuation v(x).
    # Var[v(x)] approx (epsilon^2)^2 * Var[Sum(G0(x,zi))]
    # For N i.i.d points, Var[Sum(G0)] = N * Var(G0).
    # N = 1/epsilon - 1, which we approximate as L.
    N = L
    Var_v = epsilon**4 * N * Var_G0
    
    print("Step 4: Calculate the variance of the fluctuation v(x).")
    print(f"Var[v(x)] approx epsilon^4 * N * Var[G0(x, z)], with N=L.")
    # Substitute L = 1/epsilon
    Var_v_expr = Var_v.subs(L, 1/epsilon)
    Var_v_expr = sp.simplify(Var_v_expr)
    print(f"Var[v(x)] = {Var_v_expr}\n")
    
    # Step 5: Find the maximum of Var[v(x)] over x.
    # The expression is symmetric around L/2, so the max is at x = L/2.
    x_max = L / 2
    max_Var_v = Var_v.subs(x, x_max)
    
    # Substitute L = 1/epsilon into the final expression for max variance
    max_Var_v_final = sp.simplify(max_Var_v.subs(L, 1/epsilon))

    print("Step 5: Find the maximum of Var[v(x)]. This occurs at x = L/2.")
    print(f"Max Var[v(x)] = {max_Var_v_final}\n")

    # Step 6: Calculate R = sqrt(max Var[v(x)])
    R = sp.sqrt(max_Var_v_final)
    R_coeff, R_eps_part = R.as_coeff_Mul()
    
    print("Step 6: Calculate R = (max Var[v(x)])^(1/2).")
    print("The final result for R(epsilon) is:")
    # Using 192 from sqrt(1/192) * sqrt(epsilon)
    final_coeff_denominator = int(1/R_coeff**2)
    print(f"R = (epsilon / {final_coeff_denominator})**(1/2)")
    # Also print the numeric value
    print(f"R approx {float(R_coeff)} * epsilon**(1/2)")

solve_ode_fluctuations()