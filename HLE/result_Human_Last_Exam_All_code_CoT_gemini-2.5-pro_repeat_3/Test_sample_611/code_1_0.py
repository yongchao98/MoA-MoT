import sympy

def solve_problem():
    """
    Calculates the probability P(V > 2/3) for a Beta(2,1) distributed random variable V.
    """
    # Define the symbolic variable
    x = sympy.Symbol('x')

    # Parameters for the Beta distribution
    a_param = 2
    b_param = 1

    # The PDF for Beta(2,1) is f(x) = 2x
    pdf = 2 * x

    # The value to check against
    value = sympy.Rational(2, 3)

    print("Step 1: The problem asks for the limit of P(V_n > 2/3).")
    print("Step 2: The relative area V_n converges in distribution to a random variable V.")
    print(f"Step 3: The distribution of V is known to be Beta({a_param}, {b_param}).")
    print(f"Step 4: The PDF of Beta({a_param},{b_param}) is f(x) = 2x for x in [0, 1].")
    print(f"Step 5: We calculate P(V > {value}) by integrating the PDF from {value} to 1.")
    print(f"\nP(V > {value}) = Integral from {value} to 1 of (2x) dx")

    # Calculate the antiderivative
    antiderivative = sympy.integrate(pdf, x)
    print(f"\nThe antiderivative of 2x is {antiderivative}.")

    # Evaluate the integral
    result = sympy.integrate(pdf, (x, value, 1))
    
    # Show the calculation steps with numbers
    upper_bound = 1
    lower_bound_num = 2
    lower_bound_den = 3
    
    print(f"\nEvaluating at the limits: [{antiderivative}] from {lower_bound_num}/{lower_bound_den} to {upper_bound}")
    print("This leads to the final equation:")
    
    # Print out each number in the final equation
    print(f"({upper_bound})^2 - ({lower_bound_num}/{lower_bound_den})^2 = {upper_bound*upper_bound} - {lower_bound_num**2}/{lower_bound_den**2} = {lower_bound_den**2 - lower_bound_num**2}/{lower_bound_den**2}")
    
    print(f"\nThe final probability is {result}.")

solve_problem()