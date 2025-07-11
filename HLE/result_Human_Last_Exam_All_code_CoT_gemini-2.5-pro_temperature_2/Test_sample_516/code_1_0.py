import sympy as sp

def solve_expected_value():
    """
    This function calculates the expected value of T symbolically using sympy.
    The problem is to find E[T] where T = sum_{i=1 to m} Y_i / 3^i,
    and m is the smallest integer such that Y_m > Y_{m+1}.
    Y_i are i.i.d. U(0,1).
    """

    # Define symbols for the calculation
    i = sp.Symbol('i', integer=True, positive=True)
    x = sp.Symbol('x', real=True)

    # Step 1: Express E[T] as a sum
    # E[T] = E[sum_{i=1 to inf} (Y_i / 3^i) * I(i <= m)]
    #      = sum_{i=1 to inf} (1 / 3^i) * E[Y_i * I(i <= m)]
    # The condition i <= m is equivalent to Y_1 <= Y_2 <= ... <= Y_i.

    print("Step 1: Express E[T] as an infinite sum.")
    print("E[T] = Sum_{i=1 to infinity} (1/3^i) * E[Y_i * I(Y_1 <= ... <= Y_i)]\n")

    # Step 2: Calculate E[Y_i * I(Y_1 <= ... <= Y_i)]
    # This term is an integral of y_i over the region 0 <= y_1 <= ... <= y_i <= 1 in [0,1]^i.
    # The integral evaluates to 1 / ((i-1)! * (i+1)).
    E_Yi_indicator = 1 / (sp.factorial(i - 1) * (i + 1))
    print(f"Step 2: The expectation term E[Y_i * I(Y_1 <= ... <= Y_i)] evaluates to: {E_Yi_indicator}\n")

    # Step 3: Define the summand for E[T]
    summand = (1 / 3**i) * E_Yi_indicator
    ET_sum = sp.Sum(summand, (i, 1, sp.oo))
    print(f"Step 3: The expression for E[T] becomes the sum: {ET_sum}\n")

    # Step 4: Evaluate the sum by converting it to an integral
    # We use the identity 1/(i+1) = Integral(x^i, (x, 0, 1)).
    # This turns the sum into Integral(Sum((x/3)^i / (i-1)!), (x, 0, 1)).
    # The sum Sum_{i=1 to inf} a^i / (i-1)! = a * e^a. Here a = x/3.
    integrand = (x/3) * sp.exp(x/3)
    ET_integral = sp.Integral(integrand, (x, 0, 1))
    print(f"Step 4: The sum is evaluated by converting it to an integral: E[T] = {ET_integral}\n")

    # Step 5: Calculate the final value by evaluating the integral
    expected_value = ET_integral.doit()
    print("Step 5: Evaluating the integral gives the exact value for E[T].")
    
    # Extract numbers for the final equation output
    # The result is 3 - 2*exp(1/3)
    c1 = 3
    c2 = -2
    c3 = sp.Rational(1,3)
    
    print("\nThe final equation for the expected value of T is:")
    print(f"E[T] = {c1} + ({c2}) * exp({c3.p}/{c3.q})")
    print(f"\nSo, E[T] = {expected_value}")
    
    # For verification, print the numerical value
    print(f"\nNumerical approximation: {expected_value.evalf()}")

solve_expected_value()