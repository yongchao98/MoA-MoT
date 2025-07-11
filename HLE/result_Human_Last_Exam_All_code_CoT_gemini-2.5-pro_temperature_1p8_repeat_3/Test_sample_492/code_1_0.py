import sympy

def solve_critical_time():
    """
    This function calculates the critical time 'c' for the emergence of the
    giant component in the described dynamic random graph model.
    It uses the symbolic math library sympy to perform the derivation.
    """
    # Announce the plan
    print("Finding the critical time 'c' for the emergence of the giant component.")
    print("The giant component emerges when the average degree of the graph reaches 1.")
    print("-" * 50)

    # Define symbolic variables
    t = sympy.Symbol('t', positive=True)      # General time variable
    n = sympy.Symbol('n', positive=True)      # Rate of vertex arrival (large)
    tau = sympy.Symbol('tau', positive=True)  # Integration variable for time
    c = sympy.Symbol('c', positive=True)      # Critical time we want to find

    # Step 1: Define the number of vertices V(t) at time t
    # V(t) is approximately n*t in the large n limit.
    V_t = n * t
    print(f"1. The number of vertices at time t is V(t) = {V_t}")

    # Step 2: Define the rate of edge formation, dE/d(tau)
    # The number of pairs of vertices at time tau is V(tau)^2 / 2
    # The rate of edge formation per pair is 1/n
    V_tau = n * tau
    num_pairs = V_tau**2 / 2
    rate_per_pair = sympy.Rational(1, n)
    dE_dtau = num_pairs * rate_per_pair
    
    # Simplify the expression for the rate of edge growth
    dE_dtau_simplified = sympy.simplify(dE_dtau)
    print(f"2. The rate of edge formation at time tau is dE/dtau = {dE_dtau_simplified}")

    # Step 3: Find the number of edges E(t) by integrating the rate from 0 to t
    E_t = sympy.integrate(dE_dtau_simplified, (tau, 0, t))
    print(f"3. Integrating the rate, the number of edges at time t is E(t) = {E_t}")

    # Step 4: Calculate the average degree lambda(t) = 2 * E(t) / V(t)
    avg_degree = 2 * E_t / V_t
    avg_degree_simplified = sympy.simplify(avg_degree)
    print(f"4. The average degree at time t is lambda(t) = 2*E(t)/V(t) = {avg_degree_simplified}")

    # Step 5: Solve for the critical time c where the average degree is 1
    critical_equation = sympy.Eq(avg_degree_simplified.subs(t, c), 1)

    print("\n" + "-" * 50)
    print("5. To find the critical time 'c', we set the average degree to 1.")
    
    # Print the equation with all its numbers as requested
    # The left hand side of the equation is c**2/3
    lhs = critical_equation.lhs
    # The right hand side is 1
    rhs = critical_equation.rhs
    print(f"   The equation is: {lhs} = {rhs}")

    # Solve the equation for c
    solutions = sympy.solve(critical_equation, c)

    # Since time c must be positive, we take the positive solution
    positive_solution = solutions[0]
    
    print(f"\n   Solving for c, the exact value is: c = {positive_solution}")

if __name__ == '__main__':
    solve_critical_time()