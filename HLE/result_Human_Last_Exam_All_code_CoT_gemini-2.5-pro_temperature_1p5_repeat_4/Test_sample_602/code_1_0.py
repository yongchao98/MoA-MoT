import math

def calculate_l(n: int) -> float:
    """
    Calculates the value of the function l(n) for a given integer n >= 5.

    The formula for l(n) is derived as:
    l(n) = 2 + 2/n^2 - (4*n - 2) * sqrt(n^2 - n + 1) / n^2
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    # Matrix M(n) parameters
    # a_sq = 1 - (n - 1) / (n**2)
    # b = 1 / n
    # a = math.sqrt(a_sq)
    # The direct formula is used for simplicity and to avoid intermediate steps in the code.
    
    term_sqrt = math.sqrt(n**2 - n + 1)
    
    numerator = (4 * n - 2) * term_sqrt
    
    l_n = 2 + 2 / (n**2) - numerator / (n**2)
    
    return l_n

def solve():
    """
    Solves the problem by calculating l(n) for n=5, the smallest valid integer.
    The final expression for l(n) is:
    l(n) = 2 + 2/n^2 - ( (4*n - 2) * sqrt(n*n - n + 1) ) / (n*n)
    We will print the values of the components of this equation for n=5.
    """
    n = 5
    
    # Calculate each number in the final equation for n=5
    val_2_1 = 2
    val_2_2 = 2
    val_4 = 4
    val_2_3 = -2
    val_1_1 = 1
    val_1_2 = -1
    val_1_3 = 1
    
    l_5 = calculate_l(n)

    print(f"For n = {n}, the exact expression is:")
    print(f"l({n}) = {val_2_1} + {val_2_2}/{n}^2 - (({val_4}*n + ({val_2_3})) * sqrt({val_1_1}*n^2 + ({val_1_2})*n + {val_1_3})) / {n}^2")
    print(f"l({n}) = 2 + 2/{25} - (18 * sqrt(21)) / 25")
    print(f"l({n}) = (52 - 18 * sqrt(21)) / 25")
    print(f"The numerical value is: l({n}) = {l_5}")

solve()