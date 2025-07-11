import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def get_partitions(n):
    """
    Generates all integer partitions of a non-negative integer n.
    For example, for n=3, it yields (3,), (2, 1), (1, 1, 1).
    """
    partitions = set()
    # A recursive helper function to find partitions
    def find_partitions(current_sum, start_num, terms):
        if current_sum == n:
            partitions.add(tuple(sorted(terms, reverse=True)))
            return
        if current_sum > n:
            return

        # Iterate from the last used number to n to build the partition
        for i in range(start_num, n - current_sum + 1):
            find_partitions(current_sum + i, i, terms + (i,))

    find_partitions(0, 1, tuple())
    return sorted(list(partitions), reverse=True)


def find_nonabelian_filled_groups(q, m):
    """
    Finds and describes the nonabelian filled groups of order 2*q^m.

    Args:
        q (int): An odd prime number.
        m (int): A natural number (positive integer).
    """
    # 1. Input validation
    if not isinstance(q, int) or not isinstance(m, int):
        print("Error: q and m must be integers.")
        return
    if m < 1:
        print("Error: m must be a natural number (m >= 1).")
        return
    if q % 2 == 0 or not is_prime(q):
        print(f"Error: q must be an odd prime number. {q} is not.")
        return

    order = 2 * (q**m)
    print(f"Finding nonabelian filled groups of order {order} = 2 * {q}^{m}:\n")

    # 2. Handle the known exception from group theory
    if q == 3 and m == 1:
        print("There are no nonabelian filled groups of order 6.")
        print("The group S_3 (isomorphic to Dih(Z_3)) has this structure but is the one exception and is not filled.")
        return

    # 3. Find all abelian groups P of order q^m
    # This is equivalent to finding all integer partitions of m
    partitions = get_partitions(m)
    
    if not partitions:
         print(f"Could not compute partitions for m={m}.")
         return

    print(f"The nonabelian filled groups of order {order} are the generalized dihedral groups Dih(P), where P is an abelian group of order {q**m}.")
    print("The possible structures for P are determined by the partitions of m.")
    print(f"Found {len(partitions)} such group structure(s):\n")

    # 4. Describe each group based on the structure of P
    for i, p_partition in enumerate(partitions):
        p_components = []
        for part in p_partition:
            if part == 1:
                p_components.append(f"Z_{q}")
            else:
                p_components.append(f"Z_{{{q}^{part}}}")
        
        p_structure_str = " x ".join(p_components)
        
        # Output each number in the final equation as requested
        # e.g., print "Dih(Z_{3^2})" instead of "Dih(Z_9)" to show the original numbers m_i
        group_name = f"Dih({p_structure_str})"
        
        print(f"{i+1}. {group_name}")

if __name__ == '__main__':
    # --- Examples ---

    # Case from the problem statement: q=3, m=2. Order = 2 * 3^2 = 18
    find_nonabelian_filled_groups(q=3, m=2)
    print("-" * 40)
    
    # The exceptional case: q=3, m=1. Order = 6.
    find_nonabelian_filled_groups(q=3, m=1)
    print("-" * 40)
    
    # Another example: q=5, m=3. Order = 2 * 5^3 = 250
    find_nonabelian_filled_groups(q=5, m=3)
    print("-" * 40)
    
    # An invalid input example
    find_nonabelian_filled_groups(q=9, m=2)
