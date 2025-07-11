import sympy

def find_critical_time():
    """
    This script calculates the exact time 'c' for the emergence of the
    giant component in a dynamic random graph model.

    The model is defined as:
    1. Vertices appear at a stochastic rate of n.
    2. Any possible edge between existing vertices appears at a stochastic rate of 1/n.

    The calculation is performed in the n -> infinity limit.
    """
    # Define symbolic variables for time (t) and the rate multiplier (n)
    t, n = sympy.symbols('t n', positive=True)

    # --- Step 1: Derive the average degree k(t) ---
    print("Step 1: Determine the average degree k(t) as a function of time t.")

    # In the large n limit, we can use expected values.
    # Expected number of vertices V(t) at time t, arriving at rate n.
    V_t = n * t
    print(f"Expected number of vertices V(t) = n*t")

    # The rate of change of the expected number of edges E(t) is:
    # dE/dt = (rate per edge) * (expected number of pairs of vertices)
    # dE/dt = (1/n) * E[V(t) choose 2]
    # For a Poisson process, E[V(t)*(V(t)-1)] = E[V(t)]^2 = (n*t)^2.
    # So, dE/dt = (1/n) * (n*t)**2 / 2 = n * t**2 / 2.
    
    # Integrate the rate to get the expected number of edges E(t).
    E_t = sympy.integrate(n * t**2 / 2, (t, 0, t))
    print(f"Expected number of edges E(t) = {E_t}")

    # The average degree k(t) = 2 * E(t) / V(t)
    k_t = 2 * E_t / V_t
    k_t_simplified = sympy.simplify(k_t)
    print(f"Average degree k(t) = 2*E(t)/V(t) = {k_t_simplified}")
    print("-" * 40)

    # --- Step 2: Set average degree to 1 to find critical time c ---
    print("Step 2: The giant component emerges when the average degree equals 1.")
    
    # We are looking for the critical time c where k(c) = 1.
    c = sympy.Symbol('c', positive=True)
    critical_equation = sympy.Eq(k_t_simplified.subs(t, c), 1)

    # The problem asks to output each number in the final equation.
    # The equation is c**2 / 3 = 1.
    lhs = critical_equation.lhs
    rhs = critical_equation.rhs
    base = lhs.args[0].args[0]
    exponent = lhs.args[0].args[1]
    denominator = 1/lhs.args[1]
    
    print(f"The critical equation is: {base}^{exponent} / {denominator} = {rhs}")
    print("-" * 40)

    # --- Step 3: Solve the equation for c ---
    print("Step 3: Solving the equation for the critical time c.")

    # Multiply both sides by the denominator
    step_2_equation = sympy.Eq(base**exponent, rhs * denominator)
    print(f"Multiplying by {denominator} gives: {base}^{exponent} = {rhs * denominator}")
    
    # Take the square root
    print(f"Taking the square root gives: {base} = sqrt({rhs * denominator})")

    # Use sympy's solver
    solutions = sympy.solve(critical_equation, c)
    # Since we defined c as positive, we get the single positive solution
    final_c_value = solutions[0]
    
    print("-" * 40)
    print(f"The exact value for the critical time c is: {final_c_value}")
    print(f"The approximate numerical value is: {final_c_value.evalf()}")

if __name__ == '__main__':
    find_critical_time()