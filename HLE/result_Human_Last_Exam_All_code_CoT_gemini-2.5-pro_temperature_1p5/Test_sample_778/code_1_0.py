import math

def solve_groupoid_cardinality():
    """
    Calculates the groupoid cardinality of the inserter Ins(F,G) for the given parameters.

    The groupoid cardinality of the inserter of two functors between deloopings of finite abelian groups
    is given by the ratio of the order of the codomain group to the order of the domain group.
    """
    n = 54669191328000
    m = 1470432000
    
    # The values of a and b are not needed for the cardinality calculation,
    # as the result simplifies to m/n. We can verify that the homomorphisms are well-defined:
    # a*n must be divisible by m, and b*n must be divisible by m.
    # n/m = 54669191328000 / 1470432000 = 37180.
    # So n is a multiple of m, which means a*n and b*n are definitely divisible by m.

    # The groupoid cardinality is m/n. We simplify this fraction.
    common_divisor = math.gcd(m, n)
    
    numerator = m // common_divisor
    denominator = n // common_divisor
    
    print("The groupoid cardinality of the inserter Ins(F,G) is calculated as the ratio m/n.")
    print("The final equation is:")
    print(f"{m} / {n} = {numerator} / {denominator}")

solve_groupoid_cardinality()