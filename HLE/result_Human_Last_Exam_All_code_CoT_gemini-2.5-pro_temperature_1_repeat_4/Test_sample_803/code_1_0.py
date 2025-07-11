import sympy
from sympy.ntheory import factorint
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup

def get_q_m(n):
    """Checks if n = 2 * q^m for an odd prime q. Returns (q, m) or (None, None)."""
    if n % 2 != 0 or n < 2:
        return None, None
    k = n // 2
    if k == 1:  # Group is Z_2, which is abelian
        return None, None
    factors = factorint(k)
    if len(factors) != 1:
        return None, None
    q = list(factors.keys())[0]
    m = factors[q]
    if q == 2: # q must be an odd prime
        return None, None
    return q, m

def check_if_filled_2q_m(G, group_name):
    """
    Checks if a group G is a nonabelian filled group of order 2*q^m.
    Prints a detailed analysis.
    """
    print(f"--- Checking group: {group_name} ---")
    n = G.order()
    print(f"Order = {n}")

    q, m = get_q_m(n)
    if q is None:
        print(f"Result: The group order is not of the form 2*q^m for an odd prime q.")
        print("-" * (len(group_name) + 21))
        return

    print(f"The order is of the form 2 * q^m with q = {q}, m = {m}. The equation is {n} = 2 * {q}^{m}.")

    if G.is_abelian():
        print(f"Result: The group is abelian, so it is not a NONABELIAN filled group of this type.")
        print("-" * (len(group_name) + 21))
        return
    
    print("The group is nonabelian.")
    
    # The condition for being a filled group of this type is G' = Q,
    # where G' is the commutator subgroup and Q is the Sylow q-subgroup.
    G_prime = G.commutator()
    sylow_q = G.sylow_subgroup(q)

    print(f"Order of commutator subgroup G': |G'| = {G_prime.order()}")
    print(f"Order of Sylow q-subgroup Q: |Q| = {sylow_q.order()}")

    # Since G' is a normal subgroup of G and its order is a power of q, G' must be a subgroup of Q.
    # Therefore, G' = Q if and only if their orders are equal.
    is_filled = (G_prime.order() == sylow_q.order())

    if is_filled:
        print("The condition |G'| = |Q| is met.")
        print(f"Result: {group_name} IS a nonabelian filled group of order 2*q^m.")
    else:
        print("The condition |G'| = |Q| is NOT met.")
        print(f"Result: {group_name} IS NOT a nonabelian filled group of order 2*q^m.")
    print("-" * (len(group_name) + 21))
    print()


if __name__ == '__main__':
    # Test Case 1: Dihedral group D_5 (order 10 = 2*5^1)
    D5 = sympy.DihedralGroup(5)
    check_if_filled_2q_m(D5, "Dihedral group D_5")

    # Test Case 2: Dihedral group D_9 (order 18 = 2*3^2)
    D9 = sympy.DihedralGroup(9)
    check_if_filled_2q_m(D9, "Dihedral group D_9")
    
    # Test Case 3: Dih(Z_3 x Z_3) (order 18 = 2*3^2)
    # This is G = <a,b,c | a^3=b^3=c^2=1, ab=ba, cac=a^-1, cbc=b^-1>
    a = Permutation(0, 1, 2)
    b = Permutation(3, 4, 5)
    c = Permutation(1, 2)(4, 5)
    Dih_Z3xZ3 = PermutationGroup(a, b, c)
    check_if_filled_2q_m(Dih_Z3xZ3, "Dih(Z_3 x Z_3)")

    # Test Case 4: S_3 x Z_3 (order 18 = 2*3^2), a non-filled example
    s3_gens = sympy.DihedralGroup(3).generators
    z3_gen = Permutation(3, 4, 5)
    S3xZ3 = PermutationGroup(s3_gens[0], s3_gens[1], z3_gen)
    check_if_filled_2q_m(S3xZ3, "S_3 x Z_3")

    # Test Case 5: Another non-abelian group of order 18 (SmallGroup(18,3))
    # This is G = <a,b,c | a^3=b^3=c^2=1, ab=ba, cac=a, cbc=b^-1>
    a = Permutation(0, 1, 2)
    b = Permutation(3, 4, 5)
    c = Permutation(4, 5)
    SmallGroup18_3 = PermutationGroup(a, b, c)
    check_if_filled_2q_m(SmallGroup18_3, "A non-filled group of order 18")
