import sympy

def solve_for_c():
    """
    This function calculates the critical time 'c' for the emergence of the
    giant component in the described random graph model.
    """

    # Define the symbols we will use for the calculation.
    # 'c' is the time, which must be positive.
    # 'x' is a variable for integration.
    # 'n' is the rate parameter, which tends to infinity.
    c = sympy.Symbol('c', positive=True)
    x = sympy.Symbol('x')
    n = sympy.Symbol('n', positive=True)

    print("Step 1: Define the probability distribution for the arrival of vertices.")
    print("The arrival times of any two vertices, u and v, are i.i.d. Uniform[0, c].")
    print("We need the distribution of Z = max(u, v).")
    
    # The CDF of Z is P(Z <= x) = P(u<=x, v<=x) = (x/c) * (x/c) = x**2 / c**2.
    # The PDF is the derivative of the CDF with respect to x.
    pdf_max_uv = sympy.diff((x/c)**2, x)
    print(f"The PDF of Z = max(u, v) is f(z) = {pdf_max_uv} for z in [0, c].\n")

    print("Step 2: Calculate the expected value of Z = max(u, v).")
    # E[Z] = integral from 0 to c of x * f(x) dx
    E_max_uv = sympy.integrate(x * pdf_max_uv, (x, 0, c))
    print(f"The expected value E[max(u,v)] is: {E_max_uv}\n")

    print("Step 3: Calculate the average probability of an edge.")
    # For large n, the edge probability p between two vertices is approximately (c - max(u,v))/n.
    # The average probability is E[p] = E[(c - max(u,v))/n].
    avg_p = (c - E_max_uv) / n
    print(f"The average edge probability p_avg is approximately: {avg_p}\n")

    print("Step 4: Calculate the average degree of the graph at time c.")
    # The number of vertices at time c is N ≈ n*c.
    # The average degree λ(c) ≈ N * p_avg.
    N = n * c
    avg_degree = sympy.simplify(N * avg_p)
    print(f"The average degree λ(c) is approximately: {avg_degree}\n")

    print("Step 5: Solve for the critical time 'c' where the average degree is 1.")
    # The giant component emerges when the average degree equals 1.
    critical_equation = sympy.Eq(avg_degree, 1)
    print(f"The critical equation is λ(c) = 1, which is: {critical_equation}")
    
    # We output the components of the final equation to be explicit.
    numer, denom = avg_degree.as_numer_denom()
    print(f"In the equation {numer}/{denom} = 1, the numerator is {numer}, the denominator is {denom}, and the right hand side is 1.")

    # Solve the equation for c.
    solution = sympy.solve(critical_equation, c)
    # The solution is a list, we take the single positive result.
    final_c = solution[0]
    print(f"\nSolving for c, the exact value is: {final_c}")

solve_for_c()