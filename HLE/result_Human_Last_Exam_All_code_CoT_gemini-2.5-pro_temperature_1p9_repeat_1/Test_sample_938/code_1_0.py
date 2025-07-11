import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim(ln(s_n)/n).
    The limit is determined by the dominant root of the characteristic equation
    of the recurrence relation for a_n.
    """
    
    # The value of the dominant root, lambda = 2 + sqrt(7).
    dominant_root = 2 + math.sqrt(7)
    
    # The limit is 2 * ln(lambda).
    limit_value = 2 * math.log(dominant_root)
    
    # The quantity we need to evaluate.
    final_value = 10000 * limit_value
    
    # The problem asks to output the numbers in the final equation.
    # The final equation is 10000 * 2 * ln(2 + sqrt(7)).
    val_10000 = 10000
    val_2 = 2
    sqrt_7 = math.sqrt(7)
    
    print(f"The equation to be calculated is: {val_10000} * {val_2} * ln(2 + {sqrt_7})")
    print(f"The result of the calculation is: {final_value}")
    
    integer_part = int(final_value)
    print(f"The integer part of the result is: {integer_part}")

solve()