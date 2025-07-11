import math

def count_manifolds():
    """
    Calculates the number of closed orientable 3-manifolds with
    fundamental group of order 10!.

    The solution relies on the following facts from topology and group theory:
    1. A closed orientable 3-manifold with a finite fundamental group G must be
       homeomorphic to a spherical space form S^3/G.
    2. A finite group G can act freely on S^3 if and only if its Sylow p-subgroups
       for odd primes p are cyclic, and its Sylow 2-subgroup is either cyclic
       or generalized quaternion.
    3. The order of the group is 10!, which has the prime factorization 2^8 * 3^4 * 5^2 * 7^1.
    4. Any group G satisfying these conditions for order 10! can be shown to be a direct
       product G = P_2 x H, where P_2 is the Sylow 2-subgroup and H is the subgroup of
       odd order |H| = 3^4 * 5^2 * 7^1.

    The total number of manifolds is the product of the number of choices for P_2 and H.
    """

    # Step 1: Count the number of possibilities for the Sylow 2-subgroup (P_2).
    # The order of P_2 is 2^8 = 256.
    # The condition states P_2 must be cyclic (Z_256) or generalized quaternion (Q_256).
    num_p2_choices = 2
    print(f"The Sylow 2-subgroup (order 256) can be cyclic or generalized quaternion.")
    print(f"Number of choices for the Sylow 2-subgroup = {num_p2_choices}\n")


    # Step 2: Count the number of possibilities for the odd-order part (H).
    # H is a group of order 3^4 * 5^2 * 7^1 with cyclic Sylow subgroups.
    # A group with all Sylow subgroups cyclic is metacyclic.
    # Using Sylow theorems, we can show that such a group H must be a semidirect product
    # of the form (Z_81 x Z_25) ⋊ Z_7, which simplifies to Z_2025 ⋊ Z_7.
    # The number of non-isomorphic semidirect products is determined by the number of
    # conjugacy classes of homomorphisms from Z_2025 to Aut(Z_7).
    # |Aut(Z_7)| = phi(7) = 6.
    # The order of the image of a homomorphism from Z_2025 to Z_6 must divide gcd(2025, 6) = 3.
    # So the image can be trivial (order 1) or a group of order 3.
    # This gives two non-isomorphic group structures for H.
    num_h_choices = 2
    print(f"The odd-order part H (order 14175) must have cyclic Sylow subgroups.")
    print(f"The number of non-isomorphic groups H satisfying this is 2.\n")


    # Step 3: Calculate the total number of manifolds.
    # The total number is the product of the number of choices for P_2 and H.
    total_manifolds = num_p2_choices * num_h_choices

    print(f"The total number of manifolds is the product of these choices.")
    # The final print statement shows each number in the equation as requested.
    print(f"Final calculation: {num_p2_choices} * {num_h_choices} = {total_manifolds}")
    
    return total_manifolds

if __name__ == '__main__':
    count_manifolds()