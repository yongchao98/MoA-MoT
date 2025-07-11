import sympy

# This script calculates the constant 'c', the time of emergence of the
# giant component in a variant Erdos-Renyi random graph model.
# The derivation is performed using the symbolic math library sympy.

def solve_for_c():
    """
    Derives the value of the constant c step-by-step.
    """
    # Define symbolic variables for time, n, and the constant c
    t, n, c = sympy.symbols('t n c', positive=True, real=True)
    tau = sympy.symbols('tau', positive=True, real=True)

    print("Step-by-step derivation of the constant c:")
    print("="*50)

    # Step 1: Model the expected number of vertices v(t).
    # Vertices appear at a rate 'n', so V(t) follows a Poisson process with rate n.
    # The expected number of vertices at time t is E[V(t)] = n*t.
    v_t = n * t
    print(f"1. The expected number of vertices at time t is v(t) = E[V(t)] = {v_t}.")

    # Step 2: Model the rate of change of the expected number of edges, de(t)/dt.
    # The rate is E[C(V(t), 2)] / n = E[V(t)*(V(t)-1) / (2*n)].
    # For a Poisson distribution V(t) with mean lambda = n*t, we have E[V(t)] = lambda
    # and Var(V(t)) = lambda.
    # Since Var(V(t)) = E[V(t)^2] - (E[V(t)])^2, we get E[V(t)^2] = Var + E^2.
    E_V_t_squared = v_t + v_t**2
    de_dt_expr = (E_V_t_squared - v_t) / (2 * n)
    de_dt = sympy.simplify(de_dt_expr)
    print(f"\n2. The rate of change of expected edges is de/dt = {de_dt}.")

    # Step 3: Integrate de/dt to find the expected number of edges e(t).
    e_t = sympy.integrate(de_dt.subs(t, tau), (tau, 0, t))
    print(f"\n3. Integrating de/dt gives the expected number of edges e(t) = {e_t}.")

    # Step 4: Calculate the average degree k(t).
    # k(t) is approximated by 2 * e(t) / v(t) for large n.
    k_t = sympy.simplify(2 * e_t / v_t)
    print(f"\n4. The average degree at time t is k(t) = 2 * e(t) / v(t) = {k_t}.")

    # Step 5: Solve for the critical time c.
    # The giant component emerges when the average degree k(c) equals 1.
    equation = sympy.Eq(k_t.subs(t, c), 1)
    print(f"\n5. The emergence condition is k(c) = 1, giving the equation: {equation}.")

    solutions = sympy.solve(equation, c)
    # The solution must be positive, as time is positive.
    final_c = solutions[0]

    print(f"\n6. Solving for c, we find the exact value c = {final_c}.")
    print(f"   The numerical value is c ≈ {final_c.evalf()}.")
    print("="*50)

    # As requested, outputting each number in the final equation.
    # The equation is: c^2 / 3 = 1
    power = 2
    denominator = 3
    rhs = 1
    print("\nThe final equation is c^2 / 3 = 1.")
    print("The numbers in this equation are:")
    print(f"The power of c is: {power}")
    print(f"The denominator is: {denominator}")
    print(f"The right-hand side is: {rhs}")


if __name__ == '__main__':
    solve_for_c()
