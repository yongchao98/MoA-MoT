import sympy

def solve_for_giant_component_time():
    """
    Calculates the time 'c' for the emergence of the giant component
    in the described variant of the Erdos-Renyi random graph model.
    """
    # Define symbols for our variables. We assume n and c are positive.
    n, c, s = sympy.symbols('n c s', positive=True)

    # --- Step 1: Define Expected Vertices ---
    # Vertices arrive via a Poisson process with a total rate of n.
    # The expected number of vertices at time c, E[V(c)], is the rate multiplied by time.
    E_Vc = n * c

    print("Step 1: Calculate the expected number of vertices at time c.")
    print("Vertices arrive at rate n. The expected number of vertices E[V(c)] is:")
    print(f"E[V(c)] = {E_Vc}")
    print("-" * 50)

    # --- Step 2: Define Expected Edges ---
    # At any time s <= c, the number of vertices V(s) is a Poisson random variable
    # with mean lambda = n*s.
    # The expected number of pairs of vertices is E[V(s) choose 2].
    # For a Poisson variable X with mean lam, E[X^2] = lam + lam^2.
    # E[V(s) choose 2] = E[V(s)*(V(s)-1)/2] = (E[V(s)^2] - E[V(s)])/2
    #                    = ( (n*s + (n*s)**2) - n*s ) / 2 = (n*s)**2 / 2
    E_Vs_choose_2 = (n * s)**2 / 2

    # The rate of edge formation is 1/n per pair. So, the expected rate of
    # new edge creation at time s is:
    expected_edge_rate_at_s = E_Vs_choose_2 / n

    # The total expected number of edges E[M(c)] at time c is the integral
    # of this expected rate from 0 to c.
    E_Mc = sympy.integrate(expected_edge_rate_at_s, (s, 0, c))

    print("Step 2: Calculate the expected number of edges at time c.")
    print("The expected rate of edge formation at time s is E[V(s) choose 2] / n.")
    print(f"This evaluates to: {expected_edge_rate_at_s}")
    print("Integrating this rate from 0 to c gives the expected number of edges E[M(c)]:")
    print(f"E[M(c)] = Integral({expected_edge_rate_at_s}, (s, 0, c)) = {E_Mc}")
    print("-" * 50)


    # --- Step 3: Set up and Solve the Equation ---
    # The giant component emerges when the average degree is 1, meaning 2*M/V = 1, or M = V/2.
    # We apply this to the expectations: E[M(c)] = E[V(c)] / 2.
    equation = sympy.Eq(E_Mc, E_Vc / 2)

    print("Step 3: Set up and solve the equation for c.")
    print("The giant component condition is E[M(c)] = E[V(c)] / 2.")
    print("The resulting equation is:")
    # To satisfy the prompt to show the numbers in the equation:
    # We manually format the equation string to show the expression clearly.
    print(f"{E_Mc.args[0]} * c**3 / {E_Mc.args[1]} = ({E_Vc.args[0]} * c) / 2")
    
    # Solve for c. sympy.solve returns a list of solutions.
    # Since we defined c as positive, it will only find positive solutions.
    solution = sympy.solve(equation, c)
    final_c = solution[0]

    print("\nSolving for c:")
    print(f"c**3 / 6 = c / 2")
    print(f"c**2 / 6 = 1 / 2   (since c > 0)")
    print(f"c**2 = 3")
    print(f"c = {final_c}")
    print("\n--- Final Answer ---")
    print(f"The exact value of c is the square root of 3.")
    print(f"c = {final_c.evalf()}")

if __name__ == '__main__':
    solve_for_giant_component_time()