import math

def solve_inserter_cardinality():
    """
    Calculates the groupoid cardinality of the inserter Ins(F,G) for the given problem.

    The groupoid cardinality of the inserter of two functors F, G between
    delooping groupoids of cyclic groups Z/nZ and Z/mZ, associated with
    homomorphisms f(x)=ax and g(x)=bx, is given by the ratio m/n, provided
    the homomorphisms are well-defined.

    The maps f and g are well-defined homomorphisms if n*a is divisible by m
    and n*b is divisible by m.
    """
    n = 54669191328000
    m = 1470432000
    a = 991
    b = 223
    
    # Check if the maps are well-defined homomorphisms.
    # In this problem, n is a multiple of m, so n*k % m is always 0 for any integer k.
    is_well_defined = (n % m == 0)
    
    if not is_well_defined:
        # This case is not expected based on the problem statement.
        # But as a good practice, it is checked.
        # The condition is that an and bn must be 0 mod m.
        if (n * a) % m == 0 and (n * b) % m == 0:
            is_well_defined = True

    if not is_well_defined:
        print("Error: The maps are not well-defined homomorphisms, so the problem setup is invalid.")
        return

    # The groupoid cardinality is the fraction m / n.
    # We simplify this fraction by finding the greatest common divisor (GCD).
    common_divisor = math.gcd(m, n)
    
    simplified_m = m // common_divisor
    simplified_n = n // common_divisor

    # Print the equation as requested.
    print("The groupoid cardinality is given by the formula m/n.")
    print(f"The final equation is: {m} / {n} = {simplified_m} / {simplified_n}")
    
solve_inserter_cardinality()