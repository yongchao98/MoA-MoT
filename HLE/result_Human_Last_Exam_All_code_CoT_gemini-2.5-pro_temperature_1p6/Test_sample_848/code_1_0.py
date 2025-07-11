import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim(N->inf) F(N)/ln(N).
    The limit is given by the sum of contributions from two families of solutions,
    corresponding to k=13 and k=5.
    """

    # For k=13, the characteristic equation is r^2 - 13r + 1 = 0
    # The larger root phi_13 governs the growth rate of solutions.
    phi_13 = (13 + math.sqrt(165)) / 2

    # For k=5, the characteristic equation is r^2 - 5r + 1 = 0
    # The larger root phi_5 governs the growth rate of solutions.
    phi_5 = (5 + math.sqrt(21)) / 2

    # The limit is L = 2/ln(phi_13) + 2/ln(phi_5)
    limit_value = 2 / math.log(phi_13) + 2 / math.log(phi_5)

    # We need to find the integer part of 10^4 * L
    result = 10**4 * limit_value
    
    # We are asked for the integer part of the result.
    integer_part = int(result)
    
    print(integer_part)

solve()