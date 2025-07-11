import sympy

def solve_pde_maximization():
    """
    This function solves the maximization problem step-by-step using symbolic math.
    """
    # --- Step 1: Define the problem and variables ---
    print("Step 1: Define the problem's mathematical entities symbolically.")
    # Define symbols for u(x), u_bar(x), u(x+1), and u_bar(x+1)
    u, u_bar, u_plus, u_bar_plus = sympy.symbols('u u_bar u_plus u_bar_plus', real=True)
    
    # Define the flux F and its derivatives
    F_expr = u * (1 - u)**2 * sympy.exp(-u_bar)
    F1 = sympy.diff(F_expr, u)
    F11 = sympy.diff(F1, u)
    
    print(f"Given F(u, u_bar) = {u * (1 - u)**2}*exp(-u_bar)")
    print(f"The partial derivative F_1 = dF/du = {F1}")
    print(f"The second partial derivative F_11 = d^2F/du^2 = {F11}")
    print("-" * 30)

    # --- Step 2: Derive the expression to be maximized ---
    print("Step 2: Derive the expression for E = (d/dt + F_1 * d/dx) * F_11.")
    print("This is the material derivative of F_11, which we can write as D_t(F_11).")
    print("Using the chain rule: E = (dF_11/du)*D_t(u) + (dF_11/du_bar)*D_t(u_bar)")
    
    # Derivations for D_t(u) and D_t(u_bar)
    # From PDE: u_t = -F_x = -(F_1*u_x + F_ubar*u_bar_x). F_ubar = -F.
    # D_t(u) = u_t + F_1*u_x = F*u_bar_x.
    # u_bar_x = u(x+1) - u(x) = u_plus - u.
    # D_t(u_bar) = u_bar_t + F_1*u_bar_x, where u_bar_t = -(F(x+1) - F(x)).
    
    F_plus_expr = u_plus * (1 - u_plus)**2 * sympy.exp(-u_bar_plus)
    u_bar_x = u_plus - u
    u_bar_t = -(F_plus_expr - F_expr)
    
    Dt_u = -sympy.diff(F_expr, u_bar) * u_bar_x
    Dt_u_bar = u_bar_t + F1 * u_bar_x
    
    F111 = sympy.diff(F11, u)
    F11_ubar = sympy.diff(F11, u_bar)
    
    E = F111 * Dt_u + F11_ubar * Dt_u_bar
    E = sympy.simplify(E)
    
    print("The general expression for E depends on the values u, u_bar at x, and u_plus, u_bar_plus at x+1.")
    print("-" * 30)

    # --- Step 3: Find the maximum of E on the boundary of the domain ---
    print("Step 3: Analyze E on the boundary of the domain [0, 1]^4.")
    print("Consider the boundary case: u = 0, u_bar = 0.")
    
    E_case = E.subs({u: 0, u_bar: 0})
    E_case = sympy.simplify(E_case)
    
    print(f"Substituting u=0, u_bar=0 into E gives: E = {E_case}")
    print("Let g(u_plus, u_bar_plus) = u_plus * (1 - (1-u_plus)**2 * exp(-u_bar_plus)). Then E = 4 * g.")
    print("-" * 30)
    
    # --- Step 4: Maximize the simplified expression g ---
    print("Step 4: Maximize g(u_plus, u_bar_plus) over the unit square [0,1]x[0,1].")
    g = u_plus * (1 - (1-u_plus)**2 * sympy.exp(-u_bar_plus))
    
    # For a fixed u_plus in (0,1), g is an increasing function of u_bar_plus because d/du_bar_plus of -exp(-u_bar_plus) is positive.
    # Therefore, the maximum must occur at the boundary u_bar_plus = 1.
    print("For any fixed u_plus in (0,1), g is maximized when u_bar_plus = 1.")
    g_at_ubar_plus_1 = g.subs({u_bar_plus: 1})
    
    print(f"Now we maximize g(u_plus, 1) = {g_at_ubar_plus_1} for u_plus in [0, 1].")
    
    # Check the boundaries for u_plus.
    g_at_u_plus_0 = g_at_ubar_plus_1.subs({u_plus: 0})
    g_at_u_plus_1 = g_at_ubar_plus_1.subs({u_plus: 1})
    
    # The derivative wrt u_plus shows no maximum in the interior (0,1).
    print(f"At the boundary u_plus = 0, g = {g_at_u_plus_0}")
    print(f"At the boundary u_plus = 1, g = {g_at_u_plus_1}")
    
    max_g = g_at_u_plus_1
    print(f"\nThe maximum value of g is {max_g}, which occurs when u_plus = 1.")
    print("-" * 30)

    # --- Step 5: Final Result ---
    print("Step 5: Calculate the final maximum value of E.")
    max_E = 4 * max_g
    
    print(f"The maximum of E is obtained by the following equation:")
    print(f"Max E = 4 * {max_g}")
    print("Final result:")
    print(max_E)

solve_pde_maximization()
<<<4>>>