import math

def solve_inserter_cardinality():
    """
    Calculates the groupoid cardinality of the inserter Ins(F,G) for the given parameters.

    The groupoid cardinality of the inserter of two functors F, G between the deloopings
    of cyclic groups Z/nZ and Z/mZ, associated with homomorphisms f(x)=ax and g(x)=bx,
    is given by the ratio m/n, provided f and g are valid group homomorphisms.
    """
    n = 54669191328000
    m = 1470432000
    a = 991
    b = 223

    # First, we check the necessary condition for F and G to be valid functors.
    # This requires f and g to be group homomorphisms.
    # The condition is that n*a must be divisible by m, and n*b must be divisible by m.
    if (n * a) % m != 0:
        print(f"Warning: f(x)={a}x is not a valid homomorphism from Z/nZ to Z/mZ.")
        return
    if (n * b) % m != 0:
        print(f"Warning: g(x)={b}x is not a valid homomorphism from Z/nZ to Z/mZ.")
        return

    # The groupoid cardinality is m / n.
    # We will simplify this fraction.
    common_divisor = math.gcd(m, n)
    
    numerator = m // common_divisor
    denominator = n // common_divisor
    
    # Per instructions, output each number in the final equation.
    print(f"{m} / {n} = {numerator} / {denominator}")

solve_inserter_cardinality()