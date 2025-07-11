import math

def solve():
    """
    Calculates the groupoid cardinality of the inserter Ins(F,G).
    The cardinality is given by the formula m/n.
    """
    n = 54669191328000
    m = 1470432000
    a = 991
    b = 223

    # The groupoid cardinality is m/n.
    # To simplify the fraction, we divide both by their greatest common divisor.
    common_divisor = math.gcd(m, n)
    
    numerator = m // common_divisor
    denominator = n // common_divisor

    print(f"The given values are:")
    print(f"n = {n}")
    print(f"m = {m}")
    print(f"a = {a}")
    print(f"b = {b}")
    print("\nThe groupoid cardinality is given by the equation:")
    print(f"Cardinality = m / n")
    print(f"Cardinality = {m} / {n}")
    print(f"After simplification, this becomes:")
    print(f"Cardinality = {numerator} / {denominator}")

solve()
