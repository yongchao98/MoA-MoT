import math

def calculate_inserter_cardinality():
    """
    Calculates the groupoid cardinality of the inserter Ins(F,G) and
    prints the result as a simplified fraction.
    """
    # Problem parameters
    n = 54669191328000
    m = 1470432000
    a = 991
    b = 223

    # The theory establishes that the groupoid cardinality is m/n,
    # provided the maps are well-defined homomorphisms.
    # This requires (a*n) % m == 0 and (b*n) % m == 0.
    # We can verify this:
    # (991 * 54669191328000) % 1470432000 is 0.
    # (223 * 54669191328000) % 1470432000 is 0.
    # The conditions hold.

    # The cardinality is the fraction m / n.
    # We simplify this fraction by dividing the numerator and denominator
    # by their greatest common divisor (GCD).
    common_divisor = math.gcd(m, n)
    
    simplified_numerator = m // common_divisor
    simplified_denominator = n // common_divisor

    # Output the final equation as requested.
    print(f"{m} / {n} = {simplified_numerator} / {simplified_denominator}")

calculate_inserter_cardinality()