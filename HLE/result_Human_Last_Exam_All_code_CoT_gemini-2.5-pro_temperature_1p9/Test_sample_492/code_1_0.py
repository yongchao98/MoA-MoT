import sympy

def find_critical_time_c():
    """
    Calculates the exact value of c, the time of emergence of the giant connected
    component in the specified variant of the Erdos-Renyi random graph model.
    The script follows a step-by-step derivation using symbolic mathematics.
    """
    # Define symbols for time, rate parameter, and the critical time c
    s, t, n, c = sympy.symbols('s t n c', real=True, positive=True)

    print("Step 1: Define the expected number of vertices and edges.")
    # E[N(t)] = n * t
    expected_num_vertices_at_t = n * t
    print(f"The expected number of vertices at time t, E[N(t)], is: {expected_num_vertices_at_t}")

    # The rate of edge formation at time s depends on N(s), the number of vertices at that time.
    # N(s) is a Poisson variable with mean (n*s).
    # The expected rate of edge formation E[Rate(s)] = E[N(s)*(N(s)-1)/2] * (1/n)
    # Using the property that for a Poisson variable X with mean L, E[X*(X-1)] = L^2:
    # E[Rate(s)] = (n*s)^2 / (2*n) = n*s^2 / 2
    expected_rate_of_edges = (n * s**2) / 2

    # E[M(t)] is the integral of the expected rate from 0 to t
    expected_num_edges_at_t = sympy.integrate(expected_rate_of_edges, (s, 0, t))
    print(f"The expected number of edges at time t, E[M(t)], is: {expected_num_edges_at_t}\n")

    print("Step 2: Define the average degree and the condition for the giant component.")
    # The average degree is 2 * E[M(t)] / E[N(t)]
    average_degree = 2 * expected_num_edges_at_t / expected_num_vertices_at_t
    print("The average degree is 2 * E[M(t)] / E[N(t)].")
    print("The giant component emerges when the average degree equals 1.\n")

    print("Step 3: Set up and solve the equation for the critical time c.")
    # Substitute t with c and set the expression for average degree to 1
    equation_lhs = average_degree.subs(t, c)
    final_equation = sympy.Eq(equation_lhs, 1)

    print("The equation to solve for c is:")
    # The numbers in the final equation are printed below.
    # The equation is: (c**2)/3 = 1
    c_squared_term, divisor = equation_lhs.as_numer_denom()
    print(f"{c_squared_term} / {divisor} = 1")
    
    # Solve for c
    solutions = sympy.solve(final_equation, c)
    # The positive solution is the answer
    c_value = solutions[0]
    
    print(f"\nSolving this gives c^2 = {divisor}, so c = sqrt({divisor}).\n")
    print("-----------------------------------------")
    print(f"The exact value of c is: {c_value}")
    print(f"The numerical value is approximately: {c_value.evalf()}")
    print("-----------------------------------------")

if __name__ == '__main__':
    find_critical_time_c()
    # Final answer format as requested by the user
    # <<<sqrt(3)>>> # This is for my internal thinking, final answer will be printed
    # The value is sympy.sqrt(3)
    c_val = sympy.sqrt(3)
    final_answer_string = f"<<<{c_val}>>>"
    # To conform to user requests, I need to get the exact answer and not a float.
    # The user examples are 'C' and '9.8', which implies the format could be string or float.
    # 'sqrt(3)' is a more precise representation.
    final_answer_string = "<<<sqrt(3)>>>" # This is better as it is exact.
    print(final_answer_string)