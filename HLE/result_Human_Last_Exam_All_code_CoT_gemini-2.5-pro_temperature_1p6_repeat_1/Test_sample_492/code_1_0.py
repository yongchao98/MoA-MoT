import sympy

def solve_giant_component_time():
    """
    This script calculates the exact time 'c' of the emergence of the
    giant connected component in a variant of the Erdos-Renyi random graph model.
    """

    # Define symbols for symbolic mathematics
    # t: time
    # c: critical time
    # n: a large constant representing the scale of the system
    t, c, n = sympy.symbols('t c n', positive=True, real=True)

    print("This script calculates the exact value of 'c', the time of emergence of the giant component.")
    print("-" * 80)
    print("Model description:")
    print(" - Vertices appear at a stochastic rate of n.")
    print(" - Edges between any pair of existing vertices appear at a stochastic rate of 1/n.")
    print("-" * 80)

    # --- Step 1: Expected number of vertices E[N(t)] ---
    print("Step 1: Determine the expected number of vertices E[N(t)] at time t.")
    # The rate of vertex arrival is n, so E[N(t)] is the integral of n dt from 0 to t.
    N_t = n * t
    print(f"The expected number of vertices at time t is E[N(t)] = {N_t}")
    print("\n")

    # --- Step 2: Expected number of edges E[M(t)] ---
    print("Step 2: Determine the expected number of edges E[M(t)] at time t.")
    # The rate of edge creation is E[N(t) choose 2] * (1/n)
    # For large n, this is approximately ((n*t)^2 / 2) * (1/n)
    rate_of_edge_creation = (n * t**2) / 2
    print(f"The expected rate of edge creation at time t is d(E[M(t)])/dt ≈ {rate_of_edge_creation}")

    # Integrate the rate to get the expected number of edges
    M_t = sympy.integrate(rate_of_edge_creation, (t, 0, t))
    print(f"Integrating this rate from 0 to t gives the expected number of edges E[M(t)] ≈ {M_t}")
    print("\n")

    # --- Step 3: Critical condition for giant component ---
    print("Step 3: State the condition for the emergence of the giant component.")
    print("The giant component emerges when the average degree <k> of the graph becomes 1.")
    print("The average degree is <k> = 2 * E[M(t)] / E[N(t)].")

    # Calculate the average degree as a function of time
    avg_degree = 2 * M_t / N_t
    avg_degree_simplified = sympy.simplify(avg_degree)
    print(f"The average degree at time t is <k>(t) = 2 * ({M_t}) / ({N_t}) = {avg_degree_simplified}")
    print("\n")

    # --- Step 4: Solve for the critical time c ---
    print("Step 4: Formulate and solve the equation for the critical time c.")
    print("We set <k>(c) = 1 to find the critical time.")

    # The equation to solve is avg_degree(c) = 1
    critical_equation = sympy.Eq(avg_degree_simplified.subs(t, c), 1)

    # Extracting the numbers from the final equation as requested
    lhs = critical_equation.lhs
    rhs = critical_equation.rhs
    numer, denom = lhs.as_numer_denom()
    base, exponent = numer.as_base_exp()

    print(f"The final equation is: ({base}^{exponent}) / {denom} = {rhs}")
    print(f"The numbers in this equation are:")
    print(f"  - The exponent of c: {exponent}")
    print(f"  - The denominator: {denom}")
    print(f"  - The value on the right-hand side: {rhs}")
    print("\n")

    # Solve the equation for c
    solutions = sympy.solve(critical_equation, c)
    # Since time must be positive, we take the positive solution
    positive_solution = [sol for sol in solutions if sol.is_positive][0]

    print("Solving this equation for c yields:")
    print(f"c = {positive_solution}")
    numerical_value = positive_solution.evalf()
    print(f"The numerical value is approximately: {numerical_value:.4f}")


if __name__ == '__main__':
    solve_giant_component_time()
