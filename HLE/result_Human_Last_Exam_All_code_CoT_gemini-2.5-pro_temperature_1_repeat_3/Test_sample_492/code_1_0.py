import sympy

def solve_for_critical_time():
    """
    Calculates the critical time 'c' for the emergence of the giant component
    in the described dynamic random graph model.
    """
    # Define symbolic variables for time (t, s), the rate parameter (n),
    # and the critical time (c).
    t, n, s, c = sympy.symbols('t n s c', positive=True)

    # Step 1: Define the expected number of vertices at time t.
    # Vertices arrive at rate n, so E[N(t)] = n*t.
    # For large n, we approximate N(t) by its expectation.
    E_N_t = n * t
    print("Step 1: The expected number of vertices at time t is E[N(t)].")
    print(f"E[N(t)] = {E_N_t}\n")

    # Step 2: Define the rate of edge formation at an intermediate time s.
    # The number of vertices at time s is approximately n*s.
    # The number of possible edges is (n*s choose 2), which for large n*s is ~(n*s)^2 / 2.
    num_pairs_s = (n * s)**2 / 2
    # Each potential edge appears at rate 1/n.
    rate_per_edge = 1 / n
    # The total rate of edge formation at time s is:
    edge_rate_s = num_pairs_s * rate_per_edge
    
    print("Step 2: The rate of new edge formation at time s is dE[M(s)]/ds.")
    print(f"dE[M(s)]/ds = (N(s) choose 2) * (1/n) ≈ {sympy.simplify(edge_rate_s)}\n")

    # Step 3: Calculate the expected number of edges at time t by integrating the rate.
    E_M_t = sympy.integrate(edge_rate_s, (s, 0, t))
    print("Step 3: The expected number of edges at time t, E[M(t)], is the integral of the rate.")
    print(f"E[M(t)] = ∫(from 0 to t) {sympy.simplify(edge_rate_s)} ds = {E_M_t}\n")

    # Step 4: Calculate the average degree k(t) of the graph.
    # The average degree is 2 * (number of edges) / (number of vertices).
    avg_degree_t = 2 * E_M_t / E_N_t
    simplified_avg_degree = sympy.simplify(avg_degree_t)
    print("Step 4: The average degree k(t) is 2 * E[M(t)] / E[N(t)].")
    print(f"k(t) = (2 * ({E_M_t})) / ({E_N_t}) = {simplified_avg_degree}\n")

    # Step 5: Solve for the critical time c when the average degree is 1.
    # We set k(c) = 1.
    avg_degree_c = simplified_avg_degree.subs(t, c)
    equation = sympy.Eq(avg_degree_c, 1)
    
    print("Step 5: The giant component emerges when the average degree k(c) = 1.")
    print("This gives us the final equation to solve for c:")
    # The final equation is c**2 / 3 = 1
    print(f"{avg_degree_c} = {equation.rhs}")

    # Solve the equation for c. We take the positive solution as time is positive.
    solutions = sympy.solve(equation, c)
    c_value = solutions[0] # sympy was told c is positive, so it returns the positive root

    print(f"\nSolving for c gives c = {c_value}.")
    print(f"The exact value of c is the square root of 3.")

if __name__ == '__main__':
    solve_for_critical_time()