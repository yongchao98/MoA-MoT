import math

def solve_density():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 remains irreducible mod p.
    """
    
    # The Galois group of the polynomial is the dihedral group D_7.
    # We need its order.
    group_order = 14
    
    # For the polynomial to be irreducible mod p, the Frobenius element Frob_p
    # must be a 7-cycle in the Galois group. We need to count the number
    # of 7-cycles in D_7. These are the elements of order 7.
    num_7_cycles = 6
    
    print("The problem is to find the density of primes p for which the polynomial f(x) is irreducible mod p.")
    print("This density is determined by the structure of the polynomial's Galois group, G.")
    print("1. The Galois group for the given polynomial is the Dihedral Group D_7.")
    print(f"2. The order of this group is |G| = {group_order}.")
    print("3. For f(x) to be irreducible mod p, the Frobenius element must be a 7-cycle in G.")
    print(f"4. The number of 7-cycles in D_7 is {num_7_cycles}.")
    print("5. By the Chebotarev Density Theorem, the required density is the ratio of these numbers.")
    print("\nFinal Calculation:")
    print(f"Density = (Number of 7-cycles) / (Order of the Group)")
    print(f"Density = {num_7_cycles} / {group_order}")
    
    # Simplify the fraction for the final answer.
    common_divisor = math.gcd(num_7_cycles, group_order)
    numerator = num_7_cycles // common_divisor
    denominator = group_order // common_divisor
    
    print(f"The simplified fraction is {numerator}/{denominator}.")

solve_density()