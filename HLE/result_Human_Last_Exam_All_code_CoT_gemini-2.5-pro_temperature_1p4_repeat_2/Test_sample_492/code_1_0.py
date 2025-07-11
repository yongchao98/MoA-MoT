import sympy

def solve_giant_component_time():
    """
    This script calculates the exact time 'c' for the emergence of the
    giant connected component in the described dynamic random graph model.
    It uses the sympy library for symbolic mathematics.
    """
    # Define the symbols we will use in our derivation.
    # t: time variable
    # n: scaling parameter (rate of vertex arrival)
    # c: the critical time we want to find
    t, n, c = sympy.symbols('t n c', positive=True, real=True)

    # Step 1: Define the expected rate of edge formation at time t.
    # The number of vertices V(t) follows a Poisson distribution with mean nt.
    # The expected number of pairs of vertices is E[V(t) choose 2] = E[V(t)*(V(t)-1)/2].
    # For a Poisson(lambda) variable, E[X*(X-1)] = lambda^2. Here lambda = n*t.
    # So, E[V(t) choose 2] = (n*t)**2 / 2.
    # The rate of edge formation for each pair is 1/n.
    expected_edge_formation_rate = ((n**2 * t**2) / 2) * (1/n)

    # Step 2: Integrate the rate from 0 to c to get the expected number of edges E[E(c)].
    expected_num_edges = sympy.integrate(expected_edge_formation_rate, (t, 0, c))
    
    # Step 3: Define the expected number of vertices E[V(c)].
    expected_num_vertices = n * c

    print("--- Derivation of the Critical Time 'c' ---")
    print(f"Expected number of vertices at time c, E[V(c)] = {expected_num_vertices}")
    print(f"Expected number of edges at time c, E[E(c)] = {expected_num_edges}")
    
    # Step 4: Calculate the average degree k(c) at time c.
    # The giant component emerges when the average degree k = 1.
    # k = 2 * E[edges] / E[vertices]
    avg_degree = sympy.simplify(2 * expected_num_edges / expected_num_vertices)
    
    print(f"\nThe average degree k(c) = 2 * E[E(c)] / E[V(c)] simplifies to: {avg_degree}")
    
    # Step 5: Set up and solve the equation k(c) = 1.
    final_equation = sympy.Eq(avg_degree, 1)
    
    print("\nThe giant component emerges when k(c) = 1.")
    print(f"This gives the final equation to solve for c: {final_equation}")
    
    solutions = sympy.solve(final_equation, c)
    # Since c must be a positive time, we take the positive solution.
    c_solution = solutions[0]
    
    print(f"\nSolving the equation for c yields: c = {c_solution}")
    
    # Step 6: Display the final equation with the solution substituted back in.
    print("\nVerification of the solution:")
    
    # Retrieving components from the symbolic expression for clarity
    base = c
    power = avg_degree.as_powers_dict()[c]
    denominator = avg_degree.as_numer_denom()[1]

    print(f"The equation is: ({base}^{power}) / {denominator} = 1")
    print(f"Substituting c = {c_solution}:")
    print(f"({c_solution})^{power} / {denominator} = 1")
    print(f"{c_solution**power} / {denominator} = 1")
    print(f"{sympy.simplify(c_solution**power / denominator)} = 1")

solve_giant_component_time()