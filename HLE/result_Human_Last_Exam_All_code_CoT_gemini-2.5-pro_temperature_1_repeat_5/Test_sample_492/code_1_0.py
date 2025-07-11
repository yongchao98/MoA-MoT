import sympy

def solve_for_critical_time():
    """
    Calculates the critical time 'c' for the emergence of the giant component
    in the specified dynamic random graph model.

    The method is as follows:
    1. Define expressions for the expected number of vertices E[N(t)] at time t.
    2. Define the rate of increase for the expected number of edges dE[M(t)]/dt.
    3. Integrate the rate to get the expected number of edges E[M(t)].
    4. Calculate the average degree d(t) = 2 * E[M(t)] / E[N(t)].
    5. Solve for the time 'c' where the average degree d(c) = 1.
    """
    # Define symbolic variables for time (t) and the rate parameter (n)
    t = sympy.Symbol('t', positive=True)
    n = sympy.Symbol('n', positive=True)
    
    # Step 1: Expected number of vertices at time t
    # Vertices arrive via a Poisson process with rate n.
    # E[N(t)] = n * t
    E_N_t = n * t
    
    # Step 2: Rate of increase of the expected number of edges
    # At any time s, the number of vertices N(s) follows a Poisson distribution
    # with mean n*s. The expected number of pairs of vertices is E[N(s) choose 2].
    # For a Poisson variable X with mean L, E[X*(X-1)] = L^2.
    # So, E[N(s) choose 2] = E[N(s)*(N(s)-1)/2] = (n*s)^2 / 2.
    # Since each pair forms an edge at rate 1/n, the total rate of edge
    # formation at time s is (1/n) * E[N(s) choose 2].
    s = sympy.Symbol('s', positive=True)
    rate_E_M_at_s = (1 / n) * ((n * s)**2 / 2)
    
    # Step 3: Expected number of edges at time t
    # We integrate the rate from time 0 to t.
    E_M_t = sympy.integrate(rate_E_M_at_s, (s, 0, t))
    
    # Step 4: Average degree of the graph at time t
    # The average degree is d(t) = 2 * E[M(t)] / E[N(t)]
    avg_degree = sympy.simplify(2 * E_M_t / E_N_t)
    
    # Step 5: Find the critical time c by solving d(c) = 1
    # Let c be the critical time t. We set the average degree to 1.
    print("The average degree of the graph at time t is found to be t**2 / 3.")
    print("The giant component emerges when the average degree equals 1.")
    print("Let c be the critical time. The final equation is:")

    power = 2
    denominator = 3
    rhs = 1
    
    print(f"c**{power} / {denominator} = {rhs}")

    # Solve the equation for t
    critical_equation = sympy.Eq(avg_degree, 1)
    solutions = sympy.solve(critical_equation, t)
    
    # The time c must be positive
    c = solutions[1]
    
    print(f"\nSolving this equation gives c = sqrt({denominator}).")
    print(f"The exact value of c is: {c}")

solve_for_critical_time()