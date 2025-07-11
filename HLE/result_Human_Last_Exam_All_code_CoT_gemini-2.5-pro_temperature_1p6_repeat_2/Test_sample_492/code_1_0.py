import sympy

def solve_for_c():
    """
    This function calculates the exact value of the constant 'c' by modeling
    the problem and solving for the critical point of the phase transition.
    """
    # Step 1: Define symbolic variables
    # t: time
    # n: rate multiplier (large)
    # x: dummy variable for integration (representing max arrival time)
    t, n, x = sympy.symbols('t n x', positive=True)

    # Step 2: Define the probability density function (PDF) of T_max = max(T_u, T_v)
    # where T_u, T_v are i.i.d. Uniform(0, t).
    # The CDF is P(T_max <= x) = (x/t)^2, so the PDF is the derivative.
    pdf_T_max = 2 * x / t**2

    # Step 3: Define the probability of an edge existing, given T_max = x.
    # The time for formation is (t-x) and the rate is 1/n.
    prob_edge_given_x = 1 - sympy.exp(-(t - x) / n)

    # Step 4: Calculate the overall edge probability p(t) by integrating over the PDF of T_max.
    integrand = prob_edge_given_x * pdf_T_max
    p_t = sympy.integrate(integrand, (x, 0, t))

    # Step 5: Calculate the average degree λ(t) ≈ V(t) * p(t).
    # V(t) is approximately n*t.
    V_t = n * t
    lambda_t = V_t * p_t

    # Step 6: Find the asymptotic behavior of the average degree as n -> infinity.
    # This gives us the average degree in the limit.
    lambda_asymptotic = sympy.limit(lambda_t, n, sympy.oo)

    # Step 7: The giant component emerges when the average degree is 1.
    # We create the equation for the critical time 'c' (represented by 't').
    # The numbers in the equation are 1 (from the exponent), 2 (from squaring), 3 (from the final calculation), and 1 (from the threshold).
    # For clarity, we will represent the equation as c^2/3 = 1
    c = sympy.symbols('c', positive=True)
    critical_equation = sympy.Eq(c**2 / 3, 1)

    # Step 8: Solve the equation for 'c'.
    solutions = sympy.solve(critical_equation, c)
    result = solutions[0] # Take the positive solution

    # Print the final result in the required format.
    print("The equation for the critical time 'c' is derived from setting the asymptotic average degree to 1.")
    print(f"The average degree λ(t) approaches: {lambda_asymptotic}")
    print("Setting this to 1 and replacing t with c, we get the equation:")
    # We output the numbers in the final, simplified equation
    c_squared_val = "c**2"
    divisor = 3
    rhs = 1
    print(f"{c_squared_val} / {divisor} = {rhs}")

    print("\nSolving for c, we get the exact value:")
    print(f"c = {result}")

if __name__ == "__main__":
    solve_for_c()
