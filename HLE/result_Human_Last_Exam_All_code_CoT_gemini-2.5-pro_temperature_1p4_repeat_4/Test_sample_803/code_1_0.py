import math

def _find_partitions_recursive(target, max_val, current_partition, all_partitions):
    """Helper function for recursively generating integer partitions."""
    if target == 0:
        all_partitions.append(list(current_partition))
        return

    # 'i' can be at most target, and at most the max_val from the previous step
    # to keep partitions in non-increasing order, avoiding duplicate structures.
    for i in range(min(max_val, target), 0, -1):
        current_partition.append(i)
        _find_partitions_recursive(target - i, i, current_partition, all_partitions)
        current_partition.pop()

def get_partitions(n):
    """
    Returns a list of all integer partitions of n. Each partition corresponds
    to an isomorphism class of an abelian group of order p**n.
    Example: get_partitions(3) returns [[3], [2, 1], [1, 1, 1]].
    """
    if not isinstance(n, int) or n < 0:
        return []
    if n == 0:
        return [[]]
    all_partitions = []
    _find_partitions_recursive(n, n, [], all_partitions)
    return all_partitions

def find_nonabelian_filled_groups(q, m):
    """
    This function explains and identifies the nonabelian filled groups
    of order 2*q**m for a given odd prime q and natural number m.
    """
    
    order = 2 * (q**m)
    print(f"Identifying nonabelian filled groups of order 2 * {q}^{m} = {order}")
    print("=" * 70)
    
    # Step 1: Theoretical foundation for "filled groups"
    print("\n[Step 1: Analysis of the 'filled group' property]")
    print("A key theorem states a finite solvable group is 'filled' if and only if it")
    print("does not have the Klein four-group (C\u2082 \u00d7 C\u2082) as a quotient.")
    print(f"All groups of order 2*q^m = {order} are solvable (by Burnside's p\u1d43q\u1d47 theorem).")
    print("A quotient group of order 4 (like C\u2082 \u00d7 C\u2082) cannot exist for such a group,")
    print(f"because 4 does not divide the group order {order} (since q is an odd prime).")
    print("\nConclusion: ALL groups of order 2*q**m (q odd) are filled.")
    
    # Step 2: Focus on nonabelian groups and their structure
    print("\n[Step 2: Identifying the NONABELIAN Groups of this Order]")
    print("The problem reduces to finding all nonabelian groups of order 2*q**m.")
    print("Any such group G must have a normal Sylow q-subgroup, Q, of order q**m.")
    print("The group G can always be described as a semidirect product: G = Q \u22ca C\u2082.")
    
    print("\nThis group G is nonabelian if either of these conditions is met:")
    print("  1. The subgroup Q is nonabelian itself (this requires m \u2265 3).")
    print("  2. The action of C\u2082 on Q in the semidirect product is non-trivial.")

    # Step 3: Illustrate with concrete examples and families
    print("\n[Step 3: Examples and Families of Nonabelian Filled Groups]")
    
    # Family 1: Q is cyclic (exists for m >= 1)
    print("\n--- Family 1: Dihedral Groups (when Q is cyclic) ---")
    print(f"If we take Q as the cyclic group C_{{{q}^{m}}}, the only non-trivial action of C\u2082 on Q is inversion.")
    print(f"This construction yields the well-known Dihedral Group D_{{2*q**m}}.")
    presentation = f"<r, s | r**({q**m}) = 1, s**2 = 1, s*r*s = r**(-1)>"
    print(f"Group Presentation: {presentation}")
    print(f"The final relation s*r*s = r**(-1) can be considered the 'final equation'.")
    print(f"The numbers in this final equation are the generator orders and relation exponents: {q**m}, 2, and -1.")

    # Family 2: Q is abelian but not cyclic (exists for m >= 2)
    if m >= 2:
        partitions = get_partitions(m)
        non_cyclic_partitions = [p for p in partitions if len(p) > 1]
        
        if non_cyclic_partitions:
            print("\n--- Family 2: Generalized Dihedral Groups (when Q is other abelian forms) ---")
            print(f"For m={m}, Q can have other abelian structures, corresponding to the partitions of {m}.")
            
            p = non_cyclic_partitions[0]
            structure = " \u00d7 ".join([f"C_{{{q}^{d}}}" for d in p])
            print(f"For example, Q could be {structure} (based on partition {p}).")
            print("One can form the 'generalized dihedral group' Dih(Q) where the C\u2082 acts by inverting every element of Q.")
            print(f"This gives another nonabelian filled group with the structure ({structure}) \u22ca C\u2082.")

    # Family 3: Q is nonabelian (exists for m >= 3)
    if m >= 3:
        print("\n--- Family 3: Products with Nonabelian Q (requires m \u2265 3) ---")
        print(f"For m={m}, nonabelian groups Q of order {q**m} exist.")
        print("For any such nonabelian group Q, the direct product G = Q \u00d7 C\u2082 is a simple example")
        print("of a nonabelian filled group of the required order.")
        print("More examples can be formed using non-trivial semidirect products Q \u22ca C\u2082.")

    print("\n" + "=" * 70)
    print("Summary: The nonabelian filled groups of order 2*q**m are precisely all the nonabelian")
    print("groups of that order. Their general structure is G = Q \u22ca C\u2082 as described above.")

if __name__ == '__main__':
    # You can change these values to explore other cases.
    # q must be an odd prime number.
    # m must be a natural number (integer >= 1).
    q_example = 3
    m_example = 4
    
    # Validate input
    is_prime = all(q_example % i != 0 for i in range(2, int(math.sqrt(q_example)) + 1))
    if not (is_prime and q_example > 1 and q_example % 2 != 0):
        print(f"Error: q must be an odd prime. Provided q={q_example} is not.")
    elif not (isinstance(m_example, int) and m_example >= 1):
        print(f"Error: m must be a natural number (>= 1). Provided m={m_example} is not.")
    else:
        find_nonabelian_filled_groups(q_example, m_example)