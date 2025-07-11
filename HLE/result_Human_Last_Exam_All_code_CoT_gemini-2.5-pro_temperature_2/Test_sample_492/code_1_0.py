import sympy

def find_critical_time():
    """
    This function calculates the exact value of c, the time of emergence 
    of the giant connected component in the specified random graph model.
    It uses symbolic mathematics to perform the derivation.
    """
    
    # 1. Define symbolic variables for time (t), the rate parameter (n),
    #    a dummy integration variable (s), and the critical time (c).
    #    We assume they are all positive real numbers.
    t, n, s, c = sympy.symbols('t n s c', positive=True, real=True)

    # 2. Express the expected number of vertices E[V(t)] at time t.
    #    Vertices arrive as a Poisson process with rate n, so E[V(t)] = n*t.
    expected_vertices = n * t

    # 3. Derive the rate of change for the expected number of edges E[E(t)].
    #    The rate of adding edges is (1/n) * (number of pairs of existing vertices).
    #    The number of vertices V(t) follows a Poisson distribution with mean mu = n*t.
    #    For a Poisson variable V, E[V*(V-1)] = mu**2.
    #    So, E[binom(V(t), 2)] = E[V(t)*(V(t)-1)/2] = (n*t)**2 / 2.
    #    d(E[E(t)])/dt = (1/n) * E[binom(V(t), 2)]
    rate_of_change_expected_edges = (1 / n) * ((n * t)**2 / 2)
    
    # Integrate this rate from 0 to t to find E[E(t)].
    # We substitute t with a dummy variable 's' for the integration.
    expected_edges = sympy.integrate(rate_of_change_expected_edges.subs(t, s), (s, 0, t))

    # 4. Calculate the average degree lambda(t) for large n.
    #    lambda(t) is approximately 2 * E[E(t)] / E[V(t)].
    avg_degree = sympy.simplify(2 * expected_edges / expected_vertices)

    # 5. The critical point for giant component emergence occurs when the average degree is 1.
    #    Set up the equation lambda(c) = 1.
    critical_equation = sympy.Eq(avg_degree.subs(t, c), 1)

    # 6. Solve the equation for c. Since c is time, we take the positive solution.
    solutions = sympy.solve(critical_equation, c)
    c_value = solutions[0] # The positive solution

    # 7. Print the derivation and the final result clearly.
    print("Derivation of the critical time 'c' for the emergence of the giant component:")
    print("-" * 70)
    print(f"Expected number of vertices E[V(t)] = {expected_vertices}")
    print(f"Rate of change of expected edges dE[E(t)]/dt = {rate_of_change_expected_edges}")
    print(f"Integrating the rate gives the expected number of edges E[E(t)] = {expected_edges}")
    print(f"The average degree is λ(t) ≈ 2 * E[E(t)] / E[V(t)] = {avg_degree}")
    print("\nThe giant component emerges when the average degree λ reaches the critical value of 1.")
    
    # Deconstruct the final equation to print its numbers as requested.
    # The equation is c**2/3 = 1
    lhs = critical_equation.lhs
    rhs = critical_equation.rhs
    
    power_val = lhs.as_powers_dict()[c]
    divisor_val = 1 / lhs.as_coeff_Mul()[0]

    print("\nSetting λ(c) = 1, we get the final equation:")
    print(f"c**{power_val} / {divisor_val} = {rhs}")
    
    print("\nIn this equation:")
    print(f" The exponent of c is: {power_val}")
    print(f" The divisor is: {divisor_val}")
    print(f" The value on the right-hand side is: {rhs}")

    print("\nSolving for c gives the exact value:")
    print(f"c = {c_value}")
    print(f"Numerically, c ≈ {c_value.evalf()}")
    print("-" * 70)

if __name__ == '__main__':
    find_critical_time()