import sympy

def find_upper_bound():
    """
    This function calculates the exact probability that the number of good and defective
    products becomes equal, which is the least upper bound for this probability.
    The method uses the properties of Polya's Urn and Martingale theory.
    """
    # Define the symbolic variable for our calculations
    x = sympy.Symbol('x')

    # Step 1 & 2: The limiting proportion M_inf follows a Beta(2, 1) distribution.
    # The PDF of Beta(2, 1) is f(x) = 2x.
    W0, B0 = 2, 1
    pdf = sympy.gamma(W0 + B0) / (sympy.gamma(W0) * sympy.gamma(B0)) * x**(W0 - 1) * (1 - x)**(B0 - 1)

    # Step 3: Calculate the probability P(M_inf < 1/2).
    # This is the integral of the PDF from 0 to 1/2.
    P_less_than_half = sympy.integrate(pdf, (x, 0, sympy.Rational(1, 2)))

    # Step 4 & 5: Relate P(T < infinity) to P(M_inf < 1/2).
    # Let p = P(T < infinity).
    # From the law of total probability:
    # P(M_inf < 1/2) = P(M_inf < 1/2 | T < inf) * p + P(M_inf < 1/2 | T = inf) * (1-p)
    #
    # We know P(M_inf < 1/2 | T = inf) = 0 because if T=inf, M_t > 1/2 for all t.
    # We know P(M_inf < 1/2 | T < inf) = 1/2 because the process restarts from a symmetric
    # state (k,k), leading to a symmetric Beta(k,k) distribution for the limit.
    P_M_less_than_half_given_T_finite = sympy.Rational(1, 2)
    P_M_less_than_half_given_T_inf = 0

    # Step 6: Solve for p = P(T < infinity).
    p = sympy.Symbol('p')
    equation = sympy.Eq(P_less_than_half, P_M_less_than_half_given_T_finite * p + P_M_less_than_half_given_T_inf * (1 - p))
    solution = sympy.solve(equation, p)
    prob_T_finite = solution[0]

    # Print the results of the calculation and the final answer.
    print("The problem is to find the upper bound for the probability (p) that the number of good and defective products will be equal at some point.")
    print("This can be solved by finding the exact probability, which is the least upper bound.")
    print("\nHere is the final calculation:")

    print(f"Let p be the probability P(T < infinity). We have the equation:")
    print(f"P(M_inf < 1/2) = P(M_inf < 1/2 | T < infinity) * p")
    print(f"where M_inf is the limiting proportion of good products, following a Beta(2,1) distribution.")
    print(f"The calculation yields:")
    print(f"  P(M_inf < 1/2) = {P_less_than_half}")
    print(f"  P(M_inf < 1/2 | T < infinity) = {P_M_less_than_half_given_T_finite}")
    print("\nThe final equation relating the numbers is:")
    print(f"{P_less_than_half} = {P_M_less_than_half_given_T_finite} * p")
    print(f"Solving for p, we find the probability is {prob_T_finite}.")

    print(f"\nThe upper bound for the probability is {prob_T_finite}.")
    
    # Final answer format
    print(f"<<<{prob_T_finite}>>>")

find_upper_bound()