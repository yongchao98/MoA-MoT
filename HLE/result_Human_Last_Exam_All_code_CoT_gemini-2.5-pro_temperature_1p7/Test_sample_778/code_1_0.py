import math

def calculate_groupoid_cardinality():
    """
    This function calculates the groupoid cardinality of the inserter Ins(F,G)
    for the given values of n, m, a, and b.
    """
    n = 54669191328000
    m = 1470432000
    a = 991
    b = 223

    # The general formula for the cardinality is |Ins(F,G)| = m * gcd(n, m/gcd(b-a, m)) / n.
    # Since n is a multiple of m in this problem, we can simplify the computation.
    # The result simplifies to G / (n/m), where G = gcd(n, m/gcd(b-a, m)).

    # Step 1: Calculate d = b - a
    d = b - a

    # Step 2: Calculate g = gcd(|d|, m)
    g = math.gcd(abs(d), m)

    # Step 3: Calculate L = m / g. In the theory, this corresponds to the order of
    # the image of the homomorphism k -> (b-a)k.
    L = m // g

    # Step 4: Calculate G = gcd(n, L). This is the numerator of our intermediate fraction.
    G = math.gcd(n, L)

    # Step 5: Calculate k = n / m. This is the denominator of our intermediate fraction.
    k = n // m

    # Step 6: The result is the fraction G / k. We simplify it by dividing
    # both numerator and denominator by their greatest common divisor.
    common_divisor = math.gcd(G, k)

    final_numerator = G // common_divisor
    final_denominator = k // common_divisor

    # Print the numerator and denominator of the final simplified fraction.
    print(f"{final_numerator} / {final_denominator}")

calculate_groupoid_cardinality()