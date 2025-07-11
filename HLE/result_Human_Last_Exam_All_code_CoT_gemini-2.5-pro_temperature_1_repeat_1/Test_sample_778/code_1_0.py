import math
from fractions import Fraction

def solve_groupoid_cardinality():
    """
    Calculates the groupoid cardinality of the inserter Ins(F,G) for the given parameters.
    """
    # Given values from the problem
    n = 54669191328000
    m = 1470432000
    a = 991
    b = 223

    print("Problem parameters:")
    print(f"n = {n}")
    print(f"m = {m}")
    print(f"a = {a}")
    print(f"b = {b}")
    print("-" * 30)

    # --- Step 1: Calculate the number of connected components, k ---
    # This is given by k = gcd(b-a, m).
    c = b - a
    k = math.gcd(c, m)
    
    print("Step 1: Calculate the number of connected components (k).")
    print(f"k = gcd(b - a, m)")
    print(f"k = gcd({b} - {a}, {m})")
    print(f"k = gcd({c}, {m})")
    print(f"k = {k}")
    print("-" * 30)

    # --- Step 2: Calculate the size of the automorphism group, |Aut| ---
    # The size is given by |Aut| = n / gcd(n, m/k).
    
    # We need to compute m/k first. This division is exact.
    m_div_k = m // k
    
    # Then compute gcd(n, m/k).
    g = math.gcd(n, m_div_k)
    
    # Finally, compute the size of the automorphism group. This division is also exact.
    aut_group_size = n // g

    print("Step 2: Calculate the size of the automorphism group (|Aut|).")
    print("|Aut| = n / gcd(n, m / k)")
    print(f"m / k = {m} / {k} = {m_div_k}")
    print(f"gcd(n, m / k) = gcd({n}, {m_div_k}) = {g}")
    print(f"|Aut| = {n} / {g} = {aut_group_size}")
    print("-" * 30)

    # --- Step 3: Calculate the groupoid cardinality ---
    # Cardinality = k / |Aut|
    cardinality = Fraction(k, aut_group_size)

    print("Step 3: Calculate the groupoid cardinality.")
    print("Cardinality = k / |Aut|")
    # Output the final equation with the computed numbers
    print(f"Cardinality = {k} / {aut_group_size}")
    print(f"The result as a simplified fraction is: {cardinality.numerator}/{cardinality.denominator}")
    print("-" * 30)

solve_groupoid_cardinality()