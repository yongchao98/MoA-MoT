import math
from fractions import Fraction

def solve_inserter_cardinality():
    """
    Calculates the groupoid cardinality of the inserter Ins(F,G) for the given parameters.
    """
    # Given values from the problem statement
    n = 54669191328000
    m = 1470432000
    a = 991
    b = 223

    # The groupoid cardinality of Ins(F,G) is given by the formula:
    # |Ins(F,G)| = (d * gcd(n, m/d)) / n
    # where d = gcd(b-a, m).

    # Step 1: Calculate c = b - a
    c = b - a
    # We use the absolute value for the gcd calculation, as gcd is positive.
    c_abs = abs(c)

    # Step 2: Calculate d = gcd(|b-a|, m)
    # This is the number of isomorphism classes.
    d = math.gcd(c_abs, m)

    # Step 3: Calculate m/d. The automorphism group size depends on this value.
    # We use integer division as d must divide m.
    m_div_d = m // d

    # Step 4: Calculate g = gcd(n, m/d)
    g = math.gcd(n, m_div_d)

    # Step 5: Calculate the final groupoid cardinality using the formula.
    # The result is (d * g) / n
    numerator = d * g
    denominator = n
    
    # Simplify the fraction for the final answer
    final_fraction = Fraction(numerator, denominator)

    # Print the explanation and step-by-step calculation
    print("The groupoid cardinality is calculated using the formula: (d * gcd(n, m/d)) / n")
    print(f"where n = {n}, m = {m}, a = {a}, b = {b}, and d = gcd(b-a, m).\n")
    print(f"1. Difference c = b - a = {b} - {a} = {c}")
    print(f"2. Number of isomorphism classes d = gcd(|c|, m) = gcd({c_abs}, {m}) = {d}")
    print(f"3. Value m/d = {m} / {d} = {m_div_d}")
    print(f"4. We calculate gcd(n, m/d) = gcd({n}, {m_div_d}) = {g}")
    
    # Print the final equation with all numbers filled in
    print("\nFinally, we substitute these values into the formula:")
    print(f"   ({d} * {g}) / {n}")
    print(f" = {numerator} / {denominator}")
    print(f" = {final_fraction.numerator}/{final_fraction.denominator}")
    
    # Print the final equation in the requested format
    print("\nThe final equation is:")
    print(f"({d} * {g}) / {n} = {final_fraction.numerator}/{final_fraction.denominator}")


solve_inserter_cardinality()