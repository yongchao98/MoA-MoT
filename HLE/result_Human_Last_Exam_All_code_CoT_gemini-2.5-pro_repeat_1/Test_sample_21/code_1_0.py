import math

def solve_density():
    """
    Calculates the natural density of primes p for which the given polynomial
    is irreducible mod p.
    """
    # The polynomial is of degree 7, and its Galois group G is PSL(2,7).
    # We set q=7 for our calculations.
    q = 7
    degree = 7

    # 1. Calculate the order of the group G = PSL(2,q).
    # The formula for the order of PSL(2,q) for odd prime q is (q^2 - 1)*q / 2.
    group_order = (q**2 - 1) * q // 2

    # 2. Calculate the number of elements of order 7.
    # We use Sylow's Theorems. The number of Sylow 7-subgroups, n_7,
    # must divide |G|/7 and be congruent to 1 mod 7.
    # |G|/7 = 168/7 = 24. Divisors of 24 are 1, 2, 3, 4, 6, 8, 12, 24.
    # Those congruent to 1 mod 7 are 1 and 8.
    # Since PSL(2,7) is a simple group, n_7 cannot be 1. Thus, n_7 = 8.
    num_sylow_subgroups = 8
    
    # Each Sylow 7-subgroup is cyclic and has (7-1) = 6 elements of order 7.
    # These subgroups intersect only at the identity.
    elements_of_order_7 = num_sylow_subgroups * (degree - 1)

    # 3. The density is the ratio of elements of order 7 to the group order.
    numerator = elements_of_order_7
    denominator = group_order
    
    # Simplify the fraction for the final answer.
    common_divisor = math.gcd(numerator, denominator)
    simplified_numerator = numerator // common_divisor
    simplified_denominator = denominator // common_divisor

    # 4. Print the result, showing the numbers in the final equation.
    print(f"The Galois group is G = PSL(2,7).")
    print(f"The order of G is: {denominator}")
    print(f"The number of elements of order 7 (7-cycles) in G is: {numerator}")
    print(f"The natural density is the ratio of these two numbers.")
    print(f"Final Equation: {numerator} / {denominator} = {simplified_numerator}/{simplified_denominator}")

solve_density()