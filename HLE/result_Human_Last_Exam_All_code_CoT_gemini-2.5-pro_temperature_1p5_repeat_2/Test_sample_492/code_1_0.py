import sympy as sp

def solve_giant_component_time():
    """
    This script symbolically derives and calculates the critical time 'c' for the
    emergence of the giant component in the given random graph model.
    """
    # Define symbolic variables for time (t, s) and vertex arrival rate (n)
    t, s, n = sp.symbols('t s n', positive=True)

    print("Step 1: Define the expected number of vertices at time t.")
    # V(t) is a Poisson random variable with mean n*t.
    expected_vertices_t = n * t
    print(f"E[V(t)] = {expected_vertices_t}\n")

    print("Step 2: Calculate the expected rate of edge formation.")
    # At any time 's', the number of vertices V(s) is Poisson with mean n*s.
    # The expected value of V(s) is E[V(s)] = n*s.
    # The variance of V(s) is Var(V(s)) = n*s.
    # We need E[V(s)*(V(s)-1)/2].
    # E[V(s)**2] = Var(V(s)) + E[V(s)]**2 = n*s + (n*s)**2
    expected_vertex_pairs_s = (n*s + (n*s)**2 - n*s) / 2

    # The rate of edge formation for each pair is 1/n.
    expected_edge_rate = (1/n) * expected_vertex_pairs_s
    simplified_edge_rate = sp.simplify(expected_edge_rate)
    print(f"The expected rate of new edges being added at time s is: dE[E(s)]/ds = {simplified_edge_rate}\n")

    print("Step 3: Calculate the expected number of edges at time t.")
    # Integrate the rate from 0 to t to get the total expected number of edges.
    expected_edges_t = sp.integrate(simplified_edge_rate, (s, 0, t))
    print(f"E[E(t)] = {expected_edges_t}\n")

    print("Step 4: Apply the giant component condition.")
    # The giant component emerges when the average degree is 1.
    # Average Degree = 2 * (Expected Edges) / (Expected Vertices)
    average_degree = sp.simplify(2 * expected_edges_t / expected_vertices_t)
    print(f"The average degree at time t is: {average_degree}\n")

    print("Step 5: Solve for the critical time c.")
    # We set the average degree to 1 to find the critical time c.
    c = sp.Symbol('c', positive=True)
    final_equation = sp.Eq(average_degree.subs(t, c), 1)

    # Print the equation components as requested
    lhs = final_equation.lhs
    rhs = final_equation.rhs
    print(f"The final equation for c is:")
    print(f"{lhs} = {rhs}")

    # Solve the equation for c
    solutions = sp.solve(final_equation, c)
    # The solution for c must be positive since time is positive.
    c_value = solutions[0]

    print(f"\nThe exact value of c is: {c_value}")

solve_giant_component_time()