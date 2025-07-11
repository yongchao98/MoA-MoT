import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim_{N->inf} F(N)/ln(N).
    The derivation shows this value is equal to floor(20000 / ln(phi)),
    where phi is the larger root of the characteristic equation r^2 - 13r + 1 = 0.
    """

    # The final calculation is based on the equation:
    # floor( C1 * (C2 / ln((C3 + sqrt(C4)) / C5)) )
    # where the constants are defined below.
    c_10000 = 10000
    c_2_numerator = 2
    c_13 = 13
    c_165 = 165
    c_2_denominator = 2

    print("The final result is the integer part of the following expression:")
    print(f"{c_10000} * ( {c_2_numerator} / ln( ( {c_13} + sqrt({c_165}) ) / {c_2_denominator} ) )")

    # Calculate phi, the golden ratio-like constant from the recurrence.
    phi = (c_13 + math.sqrt(c_165)) / c_2_denominator

    # Calculate the limit value.
    limit_value = c_2_numerator / math.log(phi)

    # Compute the final expression.
    result = c_10000 * limit_value

    # The problem asks for the integer part of this result.
    print("\nThe integer part of the result is:")
    print(int(result))

solve()
