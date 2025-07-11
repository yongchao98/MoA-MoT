import sympy

def solve_for_critical_time():
    """
    Calculates the exact value of the critical time c for the emergence of the
    giant component in the described random graph model.
    """
    # Define symbols for our variables
    t, n, c = sympy.symbols('t n c', positive=True)

    # Step 1: Expected number of vertices at time t
    V_t = n * t

    # Step 2: Expected rate of edge formation at time t
    # As derived in the explanation, the rate of change of expected edges is n*t^2 / 2.
    edge_formation_rate = n * t**2 / 2

    # Step 3: Expected number of edges E(t) by integrating the rate
    E_t = sympy.integrate(edge_formation_rate, (t, 0, t))

    # Step 4: Average degree k(t)
    k_t = 2 * E_t / V_t
    k_t_simplified = sympy.simplify(k_t)

    # Step 5: Solve for the critical time c where k(c) = 1
    # Formulate the equation k(c) = 1
    critical_equation = sympy.Eq(k_t_simplified.subs(t, c), 1)

    # Solve the equation for c
    solutions = sympy.solve(critical_equation, c)
    # Since t and c are defined as positive, we will get a single positive solution.
    c_val = solutions[0]
    
    # Print the explanation and the steps of the final calculation
    print("The emergence of the giant component occurs when the average degree k(c) equals 1.")
    print(f"The average degree as a function of time t is k(t) = {k_t_simplified}.")
    print("\nSetting t=c and k(c)=1, we get the equation for the critical time c:")
    
    # Extract numbers for printing the equation c^2 / 3 = 1
    # The equation is c**2/3 = 1
    numerator, denominator = k_t_simplified.subs(t,c).as_numer_denom()
    rhs = critical_equation.rhs
    print(f"{numerator} / {denominator} = {rhs}")
    
    # Solve for c^2
    c_squared_val = rhs * denominator
    print(f"{numerator} = {c_squared_val}")
    
    # Final value for c
    print(f"c = sqrt({c_squared_val})")
    print(f"\nThe exact value of c is: {c_val}")


solve_for_critical_time()