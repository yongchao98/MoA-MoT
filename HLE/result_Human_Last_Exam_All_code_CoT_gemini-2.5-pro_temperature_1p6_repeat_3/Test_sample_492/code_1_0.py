import sympy

def solve_for_c():
    """
    This function calculates the time 'c' of the emergence of the giant
    connected component in the specified random graph model.
    """
    # Define symbols for our symbolic calculations.
    # n is the rate parameter, t is time. We assume they are positive.
    t, n = sympy.symbols('t n', positive=True)

    # Step 1: Define the expected number of vertices at time t.
    # Vertices appear at rate n, so V(t) ≈ n*t.
    V_t = n * t
    print(f"Step 1: The expected number of vertices at time t, V(t), is modeled as:")
    print(f"V(t) = {V_t}\n")

    # Step 2: Define the rate of edge formation (dE/dt).
    # The number of pairs of vertices is V(t) choose 2 ≈ V(t)**2 / 2 for large V(t).
    # Each pair forms an edge at rate 1/n.
    # So, dE/dt = (V(t)**2 / 2) / n.
    dE_dt_rate = (V_t**2 / 2) / n
    dE_dt_simplified = sympy.simplify(dE_dt_rate)
    print(f"Step 2: The rate of edge formation, dE/dt, is (V(t)^2 / 2) / n.")
    print(f"dE/dt = (({V_t})^2 / 2) / n = {dE_dt_simplified}\n")

    # Step 3: Calculate the expected number of edges E(t) by integrating the rate.
    # We use a dummy variable 'u' for integration.
    u = sympy.symbols('u')
    E_t = sympy.integrate(dE_dt_rate.subs(t, u), (u, 0, t))
    print(f"Step 3: The expected number of edges, E(t), is the integral of the rate from 0 to t.")
    print(f"E(t) = Integrate({dE_dt_simplified.subs(t,u)}, (u, 0, t)) = {E_t}\n")

    # Step 4: Calculate the average degree k(t).
    # The average degree is k(t) = 2 * E(t) / V(t).
    k_t = 2 * E_t / V_t
    k_t_simplified = sympy.simplify(k_t)
    print(f"Step 4: The average degree, k(t), is 2 * E(t) / V(t).")
    print(f"k(t) = 2 * ({E_t}) / ({V_t}) = {k_t_simplified}\n")

    # Step 5: Solve for the emergence time c.
    # The giant component emerges when the average degree k(c) = 1.
    c = sympy.symbols('c', positive=True)
    emergence_equation = sympy.Eq(k_t_simplified.subs(t, c), 1)
    
    # We rearrange the equation to the form c^2 = constant
    lhs = emergence_equation.lhs * 3
    rhs = emergence_equation.rhs * 3
    final_equation_lhs_base = c
    final_equation_lhs_exp = sympy.degree(lhs, c)
    final_equation_rhs = rhs

    print(f"Step 5: To find c, we solve the emergence equation k(c) = 1:")
    print(f"{emergence_equation.lhs} = {emergence_equation.rhs}")
    print(f"Multiplying both sides by 3, we get the final equation:")
    print(f"{final_equation_lhs_base}^{final_equation_lhs_exp} = {final_equation_rhs}")

    # Solve for c
    solution = sympy.solve(emergence_equation, c)[0]
    print(f"\nThe exact value of c is the square root of {final_equation_rhs}, which is: {solution}")

if __name__ == '__main__':
    solve_for_c()