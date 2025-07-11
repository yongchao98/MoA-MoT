import math

def solve_polynomial_density():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 remains irreducible mod p.

    This problem is solved using the Chebotarev Density Theorem. The density is
    the ratio of the number of 7-cycles in the polynomial's Galois group (G)
    to the order of the group.

    For the given polynomial, the Galois group G is known to be the
    Affine General Linear Group AGL(1,7).
    """

    # The order of the Galois group G = AGL(1,7) is |F_7^*| * |F_7|.
    # |F_7^*| = 6, |F_7| = 7
    group_order = 6 * 7

    # The number of 7-cycles in AGL(1,7) corresponds to the affine maps
    # x -> ax + b where a=1 and b is not 0. There are 6 such choices for b.
    num_7_cycles = 6

    # The density is the ratio of the number of 7-cycles to the group order.
    numerator = num_7_cycles
    denominator = group_order

    # Simplify the fraction by dividing by the greatest common divisor.
    common_divisor = math.gcd(numerator, denominator)
    simplified_num = numerator // common_divisor
    simplified_den = denominator // common_divisor

    print("The problem asks for the natural density of primes p where the polynomial f(x) remains irreducible mod p.")
    print("This density is determined by the structure of the polynomial's Galois group, G, via the Chebotarev Density Theorem.")
    print("\n1. The Galois group for this polynomial is G = AGL(1,7).")
    print(f"2. The order of this group is {group_order}.")
    print(f"3. The number of elements in G that are 7-cycles (which correspond to irreducible factorizations mod p) is {num_7_cycles}.")
    print("\n4. The density is the ratio of these two numbers.")
    print(f"Density = {numerator} / {denominator}")
    print(f"Simplified Density = {simplified_num} / {simplified_den}")


if __name__ == "__main__":
    solve_polynomial_density()