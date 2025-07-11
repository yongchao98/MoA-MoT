import sympy

def solve_for_c():
    """
    This script derives the exact value of 'c', the time of emergence of the
    giant component in the described variant of the Erdos-Renyi model.
    """
    print("This program will derive the value of 'c' step by step.\n")

    # Define symbolic variables for time 't' and rate multiplier 'n'.
    # 'c' will be used for the final critical time.
    t, n, c = sympy.symbols('t n c', positive=True, real=True)

    # --- Step 1: Find the expected number of vertices v(t) ---
    # Vertices appear via a Poisson process with rate 'n'.
    # The expected number of vertices at time 't' is the rate multiplied by time.
    v_t = n * t
    print("Step 1: Find the expected number of vertices, v(t).")
    print(f"         Vertices appear at rate n, so v(t) = E[V(t)] = n * t.\n")

    # --- Step 2: Set up and solve for the expected number of edges e(t) ---
    # The rate of change of the expected number of edges (de/dt) is the
    # rate per edge (1/n) multiplied by the expected number of pairs of vertices.
    # For a Poisson-distributed variable V(t) with mean lambda = n*t,
    # the expectation of V(t)-choose-2 is E[V(t)*(V(t)-1)/2] = lambda^2 / 2.
    de_dt = (1/n) * ((n*t)**2 / 2)
    de_dt_simplified = sympy.simplify(de_dt)
    print("Step 2: Find the expected number of edges, e(t).")
    print(f"         The rate of new edge formation is de/dt = (1/n) * E[V(t) C 2].")
    print(f"         de/dt = (1/n) * (n*t)^2 / 2 = {de_dt_simplified}")

    # To find e(t), we integrate de/dt from 0 to t.
    s = sympy.symbols('s') # dummy variable for integration
    e_t = sympy.integrate(de_dt_simplified.subs(t, s), (s, 0, t))
    print(f"         Integrating the rate gives e(t) = {e_t}.\n")

    # --- Step 3: Find the average degree of the graph, k(t) ---
    # The average degree is twice the number of edges divided by the number of vertices.
    k_t = 2 * e_t / v_t
    k_t_simplified = sympy.simplify(k_t)
    print("Step 3: Find the average degree of the graph, k(t).")
    print(f"         k(t) = 2 * e(t) / v(t) = 2 * ({e_t}) / ({v_t})")
    print(f"         After simplification, k(t) = {k_t_simplified}.\n")

    # --- Step 4: Solve for the critical time c ---
    # The giant component emerges when the average degree is 1.
    # We set k(c) = 1 and solve for c.
    print("Step 4: Solve for the critical time 'c' where the giant component emerges.")
    print("         This happens when the average degree k(c) equals 1.")
    print("         The final equation to solve is:")
    
    # We solve k(c) = 1, which is c^2 / 3 = 1.
    # This loop outputs each number in the equation per the instructions.
    equation_lhs = k_t_simplified.subs(t,c)
    print(f"         {equation_lhs} = 1")
    print(f"         c**2 = {3 * 1}")
    solution = sympy.sqrt(3)
    print(f"         c = sqrt(3)\n")

    print(f"The exact value of c is {solution}.")

if __name__ == '__main__':
    solve_for_c()

<<<sqrt(3)>>>