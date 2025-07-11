import sys

def find_partitions(n):
    """
    Generates all integer partitions of n.
    This is a helper function to find all possible structures of an
    abelian group of order q^m.
    """
    # Using a cache to store results of partitions(k, max_val)
    cache = {}

    def partitions_recursive(k, max_val):
        if (k, max_val) in cache:
            return cache[(k, max_val)]
        if k == 0:
            return [[]]
        if k < 0 or max_val == 0:
            return []
        
        # Partitions that do not include max_val
        res = partitions_recursive(k, max_val - 1)
        
        # Partitions that include at least one max_val
        partitions_with_max = partitions_recursive(k - max_val, max_val)
        for p in partitions_with_max:
            res.append(p + [max_val])
            
        cache[(k, max_val)] = res
        return res

    return partitions_recursive(n, n)

def describe_group(partition, q):
    """
    Creates a string description of an abelian group
    based on a partition and a prime q.
    """
    parts = sorted(partition, reverse=True)
    group_parts = []
    for part in parts:
        group_parts.append(f"Z_{q}^{part}")
    return " x ".join(group_parts)

def solve_filled_groups(q, m):
    """
    Identifies and describes the nonabelian filled groups of order 2*q^m.
    """
    if not isinstance(q, int) or not isinstance(m, int):
        print("Error: q and m must be integers.", file=sys.stderr)
        return
    if m <= 0:
        print(f"Error: m must be a natural number (m > 0), but got m={m}.", file=sys.stderr)
        return

    # Check if q is an odd prime
    is_prime = q > 1 and all(q % i for i in range(2, int(q**0.5) + 1))
    if not is_prime or q == 2:
        print(f"Error: q must be an odd prime, but got q={q}.", file=sys.stderr)
        return

    print(f"The nonabelian filled groups of order 2 * {q}^{m} are the generalized dihedral groups Dih(A), where A is an abelian group of order {q}^{m}.")
    print("The possible structures for A are determined by the partitions of m.")
    print("-" * 20)

    partitions = find_partitions(m)
    
    if not partitions or (len(partitions) == 1 and not partitions[0]):
        # This case happens for m=0, handled above, but as a safeguard.
        print(f"There are no abelian groups of order {q}^{m} for m={m}.")
        return

    print(f"For m = {m}, there are {len(partitions)} such non-isomorphic abelian groups, leading to {len(partitions)} families of solutions:")

    for i, p in enumerate(sorted(partitions, key=lambda x: (len(x), sorted(x)), reverse=True)):
        abelian_group_A = describe_group(p, q)
        print(f"{i+1}. The group Dih({abelian_group_A}).")
        if p == [m]:
            print(f"   This is the standard dihedral group D_{2 * q**m}.")
        else:
            print(f"   This is the semidirect product ({abelian_group_A}) x_inv Z_2, where the action is inversion.")

if __name__ == '__main__':
    # Example usage.
    # The problem asks for the general case, this demonstrates the script's output.
    # Let's use q=3 and m=4 as an example.
    
    # We must provide q and m. Since we cannot ask for input, we set them here.
    q_example = 3
    m_example = 4

    solve_filled_groups(q_example, m_example)
