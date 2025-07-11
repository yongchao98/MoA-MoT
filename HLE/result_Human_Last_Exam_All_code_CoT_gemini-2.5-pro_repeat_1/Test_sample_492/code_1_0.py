import sympy

def solve_for_c():
    """
    This function calculates the critical time 'c' for the emergence of the
    giant component in the described graph model.
    """
    # Define the symbolic variables needed for the calculation.
    t, n, c, x = sympy.symbols('t n c x', positive=True, real=True)

    print("Step 1: Determine the expected value of max(T_u, T_v) for two vertex arrival times.")
    # The arrival times T_u and T_v are independent and uniform in [0, t].
    # The probability density function (PDF) of X = max(T_u, T_v) is f(x) = 2*x/t**2.
    pdf_max_T = 2 * x / t**2
    # The expected value is the integral of x * f(x) from 0 to t.
    E_max_T = sympy.integrate(x * pdf_max_T, (x, 0, t))
    print(f"E[max(T_u, T_v)] = {E_max_T}\n")

    print("Step 2: Calculate the average edge probability p(t).")
    # For large n, the probability of an edge is approximately (t - max(T_u, T_v))/n.
    # We average this over the distribution of max(T_u, T_v).
    # p(t) = E[(t - max(T_u, T_v))/n] = (t - E[max(T_u, T_v)])/n
    p_t = (t - E_max_T) / n
    print(f"The average edge probability p(t) is approximately: {p_t}\n")

    print("Step 3: Calculate the average vertex degree lambda(t).")
    # The number of vertices V(t) at time t is approximately n*t.
    # The average degree lambda(t) = V(t) * p(t).
    lambda_t = (n * t) * p_t
    print(f"The average degree lambda(t) is: {lambda_t}\n")

    print("Step 4: Set up and solve the critical equation lambda(c) = 1.")
    # The giant component emerges when the average degree equals 1.
    critical_equation = sympy.Eq(lambda_t.subs(t, c), 1)
    print(f"The critical equation is: {critical_equation}\n")

    print("Step 5: Extract the numbers from the final equation and solve for c.")
    # The final equation is c**2 / 3 = 1.
    # We identify the numbers that make up this equation.
    lhs = critical_equation.lhs
    # The equation is of the form a*c**p / d = r
    # In our case: c**2 / 3 = 1
    # We can get the exponent by inspecting the expression.
    power = lhs.as_powers_dict()[c]
    # The denominator can also be extracted.
    denominator = 1 / lhs.as_coeff_Mul(c**power)[0]
    right_hand_side = critical_equation.rhs

    print("The numbers in the final equation c**2 / 3 = 1 are:")
    print(f"  - The exponent of c is: {power}")
    print(f"  - The denominator is: {denominator}")
    print(f"  - The value on the right-hand side is: {right_hand_side}\n")

    # Solve the equation for c. We take the positive solution as time is positive.
    solutions = sympy.solve(critical_equation, c)
    positive_solution = solutions[0] # sympy already filtered for positive c

    print(f"The exact value of c is the solution to the equation.")
    print(f"c = {positive_solution}")


if __name__ == '__main__':
    solve_for_c()