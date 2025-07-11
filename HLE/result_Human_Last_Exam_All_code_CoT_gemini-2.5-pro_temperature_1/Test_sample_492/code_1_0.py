import sympy

def solve_for_c():
    """
    This script calculates the exact value of the constant c, which represents
    the time of emergence of the giant connected component in the specified
    variant of the Erdos-Renyi random graph model.
    
    The method involves these theoretical steps, which are implemented using
    the symbolic mathematics library SymPy:
    
    1.  Model Setup: At any time `t`, the number of vertices `V` is a Poisson
        random variable with mean `n*t`. For large `n`, `V ≈ n*t`. The arrival
        times of any two vertices are independent and uniformly distributed on `[0, t]`.
    
    2.  Edge Probability: The probability `p` of an edge existing between two
        vertices depends on their arrival times. We calculate the average
        edge probability `p(t)` over all pairs of vertices present at time `t`.
    
    3.  Average Degree: The state of the graph is analogous to a `G(V, p)`
        random graph. The critical parameter for the emergence of a giant
        component is the average degree, which is `λ = V * p(t)`.
        
    4.  Phase Transition: In the limit as `n → ∞`, the giant component
        emerges when `λ` crosses the threshold of 1. We compute `lim_{n→∞} λ(t)`
        and set it equal to 1.
        
    5.  Solve for c: Solving the resulting equation for `t` gives the
        constant `c`.
    """
    print("Starting the symbolic calculation for the critical time c.")
    
    # Define the necessary symbolic variables
    t = sympy.Symbol('t', positive=True)
    z = sympy.Symbol('z')
    n = sympy.Symbol('n', positive=True)

    # For two vertices that have arrived by time t, their arrival times X and Y
    # are i.i.d. Uniform(0, t). Let Z = max(X, Y). The PDF of Z is f_Z(z) = 2*z/t**2.
    pdf_Z = 2 * z / t**2

    # The probability of an edge is p(t) = 1 - E[exp(-(t-Z)/n)], where the
    # expectation is over the distribution of Z.
    integrand = sympy.exp(-(t - z) / n) * pdf_Z
    expected_value = sympy.integrate(integrand, (z, 0, t))
    p_t = 1 - expected_value

    # The number of vertices V(t) is approximately n*t.
    # The average degree is λ ≈ V(t) * p(t) for large n.
    avg_degree = n * t * p_t

    # Calculate the limit of the average degree as n approaches infinity.
    lambda_limit = sympy.limit(avg_degree, n, sympy.oo)
    
    print("\nStep 1: The average degree of a vertex at time t in the n->infinity limit is:")
    print(f"λ(t) = {lambda_limit}")

    # The giant component emerges when the average degree λ equals 1.
    critical_equation = sympy.Eq(lambda_limit, 1)
    
    print("\nStep 2: The equation for the critical time c is λ(c) = 1. This gives the final equation:")
    # We explicitly print the numbers in the final equation as requested.
    num, den = lambda_limit.as_numer_denom()
    rhs = critical_equation.args[1]
    print(f"{num} / {den} = {rhs}")

    # Solve the equation for t to find the value of c.
    solutions = sympy.solve(critical_equation, t)
    
    # Time c must be positive.
    c = [sol for sol in solutions if sol.is_positive][0]
    
    print(f"\nStep 3: Solving the equation gives the exact value for c:")
    print(f"c = {c}")

if __name__ == '__main__':
    solve_for_c()
