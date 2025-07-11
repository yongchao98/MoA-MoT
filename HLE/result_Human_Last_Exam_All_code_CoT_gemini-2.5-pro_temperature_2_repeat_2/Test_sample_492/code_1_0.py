import sympy as sp

def solve_for_c():
    """
    This function carries out the derivation to find the constant c, which is
    the time of emergence of the giant component in the described random graph model.
    It uses the sympy library for symbolic mathematics.
    """

    # --- 1. Define Symbols ---
    # t: general time variable
    # s: integration variable for time
    # n: scaling parameter (rate of vertex arrival)
    # c: the critical time we want to find
    t, s, n, c = sp.symbols('t s n c', positive=True)

    print("Step 1: Define the expected number of vertices at time 's'.")
    # Vertices arrive via a Poisson process with rate n.
    # So, the number of vertices at time s, N(s), is a Poisson(n*s) random variable.
    # The expected value E[N(s)] is n*s.
    E_N_s = n * s
    print(f"The number of vertices N(s) follows a Poisson process with rate n.")
    print(f"Therefore, the expected number of vertices E[N(s)] = {E_N_s}\n")

    print("Step 2: Calculate the expected rate of edge formation.")
    # The instantaneous rate of edge creation is (N(s) choose 2) / n.
    # We need its expectation: E[N(s)*(N(s)-1) / (2*n)] = (E[N(s)**2] - E[N(s)]) / (2*n).
    # For a Poisson(lambda) variable, E[X] = lambda and Var(X) = lambda.
    # The second moment E[X**2] = Var(X) + E[X]**2.
    # For N(s), lambda = n*s. So E[N(s)**2] = n*s + (n*s)**2.
    E_N_s_sq = n * s + (n * s)**2
    expected_edge_rate = sp.simplify((E_N_s_sq - E_N_s) / (2 * n))
    print("The rate of edge addition at time s is (N(s) choose 2) / n.")
    print(f"The expected rate is E[(N(s)^2 - N(s)) / (2n)], which simplifies to: {expected_edge_rate}\n")

    print("Step 3: Calculate the expected total number of edges at time 't'.")
    # We integrate the expected rate from 0 to t to get E[E(t)].
    E_E_t = sp.integrate(expected_edge_rate, (s, 0, t))
    print("By integrating the edge rate from 0 to t, we get the expected number of edges E[E(t)].")
    print(f"E[E(t)] = Integral({expected_edge_rate}, (s, 0, t)) = {E_E_t}\n")

    print("Step 4: Calculate the expected average degree at time 't'.")
    # The average degree is approx. 2 * E[number of edges] / E[number of vertices].
    E_N_t = n * t
    avg_degree_t = sp.simplify(2 * E_E_t / E_N_t)
    print("The average degree is defined as 2 * (Number of Edges) / (Number of Vertices).")
    print("Using expected values, the average degree at time t is 2 * E[E(t)] / E[N(t)].")
    print(f"This gives an average degree of: {avg_degree_t}\n")

    print("Step 5: Solve for the critical time 'c'.")
    # The giant component emerges when the average degree is 1. We solve for c.
    # Set the average degree expression (with t=c) equal to 1.
    final_equation = sp.Eq(avg_degree_t.subs(t, c), 1)
    print("The giant component phase transition occurs when the average degree equals 1.")
    print("We set the average degree expression to 1 and solve for the critical time c.")

    # We manually format the equation to show all numbers.
    base, exp = final_equation.lhs.as_base_exp()
    if isinstance(final_equation.lhs, sp.Mul): # Handles c**2 / 3
        coeff, term = final_equation.lhs.as_coeff_Mul()
        base, exp = term.as_base_exp()
        print(f"The final equation is: ({base}**{exp}) / {1/coeff} = 1")
    else: # Handles simpler cases
        print(f"The final equation is: {base}**{exp} = 1")

    # Solve the equation for c.
    solutions = sp.solve(final_equation, c)
    # Since time c must be positive, we take the positive solution.
    c_value = solutions[0]

    print(f"\nSolving for c gives the exact value: c = {c_value}")

if __name__ == '__main__':
    solve_for_c()
