import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim (ln(s_n)/n)
    where the limit is 2 * ln(2 + sqrt(7)).
    """
    
    # The final expression to calculate is 10000 * 2 * ln(2 + sqrt(7)).
    # Let's break down the calculation and print each part as requested.

    # Numbers from the final equation
    multiplier = 10000
    factor1 = 2
    inner_val1 = 2
    inner_val2 = 7
    
    print(f"The equation is: floor({multiplier} * {factor1} * ln({inner_val1} + sqrt({inner_val2})))")
    
    # Step 1: Calculate the value of the dominant root lambda = 2 + sqrt(7)
    sqrt_of_7 = math.sqrt(inner_val2)
    lambda_val = inner_val1 + sqrt_of_7
    print(f"The value of the dominant root lambda is: {lambda_val}")
    
    # Step 2: Calculate the limit L = 2 * ln(lambda)
    limit_val = factor1 * math.log(lambda_val)
    print(f"The value of the limit L = 2 * ln(lambda) is: {limit_val}")
    
    # Step 3: Multiply by 10^4
    final_value = multiplier * limit_val
    print(f"The value of 10^4 * L is: {final_value}")
    
    # Step 4: Find the integer part
    integer_part = int(final_value)
    print(f"The integer part of the final value is: {integer_part}")

solve()