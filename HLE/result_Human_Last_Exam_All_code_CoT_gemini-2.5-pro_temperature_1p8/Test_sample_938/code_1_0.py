import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim_{n->inf} ln(s_n)/n.
    The step-by-step derivation shows the limit is ln(2 + sqrt(7)).
    This script calculates the numerical value and its integer part.
    """

    # The value of the limit L is ln(2 + sqrt(7))
    sqrt_7 = math.sqrt(7)
    val_2_plus_sqrt_7 = 2 + sqrt_7
    limit = math.log(val_2_plus_sqrt_7)
    
    # We need to find the integer part of 10^4 * L
    factor = 10000
    final_value = factor * limit
    integer_part = int(final_value)

    print("The expression to be calculated is the integer part of: 10000 * ln(2 + sqrt(7))")
    print("-" * 30)
    print("Step 1: Calculate the terms inside the expression.")
    print(f"The value of sqrt(7) is approximately: {sqrt_7}")
    print(f"The value of 2 + sqrt(7) is approximately: {val_2_plus_sqrt_7}")
    print(f"The value of ln(2 + sqrt(7)) is approximately: {limit}")
    print("-" * 30)
    print("Step 2: Calculate the final value.")
    print(f"The value of 10000 * ln(2 + sqrt(7)) is approximately: {final_value}")
    print("-" * 30)
    print("Step 3: Find the integer part.")
    print(f"The integer part of the final value is: {integer_part}")

solve()
