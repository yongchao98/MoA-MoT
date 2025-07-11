import math
from itertools import combinations

def get_partitions(n):
    """Generates all integer partitions of n."""
    a = [0] * (n + 1)
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

def describe_group_structure(partition, q):
    """Formats the description of an abelian group from its partition."""
    if not partition:
        return "{e}"
    
    # Count occurrences of each part
    counts = {}
    for part in partition:
        counts[part] = counts.get(part, 0) + 1
    
    # Format the string
    parts_str = []
    for part, count in sorted(counts.items()):
        term = f"C_{q}^{part}" if part > 1 else f"C_{q}"
        if count > 1:
            parts_str.append(f"({term})^{count}")
        else:
            parts_str.append(term)
    return " x ".join(parts_str)

def list_filled_groups_abelian_sylow(q, m):
    """
    Lists the structures of nonabelian filled groups of order 2*q^m
    for the case where the Sylow q-subgroup P is abelian.
    
    The general solution is that ALL nonabelian groups of order 2*q^m are filled.
    This function describes a large, well-structured subset of them.
    """
    print(f"The nonabelian filled groups of order 2*q^m are precisely all nonabelian groups of that order.")
    print(f"Below is a structural description for the cases where the Sylow q-subgroup P (of order q^{m}) is abelian.")
    print(f"These groups have the form G = A x D(B), where P = A x B, B is non-trivial, and D(B) is the generalized dihedral group over B.")
    print("-" * 60)
    print(f"Listing structures for q={q}, m={m}:")
    print("-" * 60)

    if not (isinstance(q, int) and q > 2 and all(q % i != 0 for i in range(2, int(math.sqrt(q)) + 1))):
        print("Error: q must be an odd prime.")
        return
    if not (isinstance(m, int) and m > 0):
        print("Error: m must be a natural number.")
        return

    # A set to store unique group structures to avoid duplicates
    # (e.g., from partitions with repeated parts)
    seen_structures = set()

    partitions = get_partitions(m)
    for p in partitions:
        factors = list(p)
        num_factors = len(factors)
        
        # Iterate through all non-empty subsets of factors for B
        for i in range(1, num_factors + 1):
            for subset_indices in combinations(range(num_factors), i):
                b_partition = [factors[j] for j in subset_indices]
                a_partition = [factors[j] for j in range(num_factors) if j not in subset_indices]
                
                b_partition.sort(reverse=True)
                a_partition.sort(reverse=True)

                # B must be non-trivial
                if not b_partition:
                    continue

                a_desc = describe_group_structure(a_partition, q)
                b_desc = describe_group_structure(b_partition, q)
                
                # The final equation for the group structure
                if a_desc == "{e}":
                    group_desc = f"G = D({b_desc})"
                else:
                    group_desc = f"G = {a_desc} x D({b_desc})"
                
                # Order of A and B
                a_order_m = sum(a_partition)
                b_order_m = sum(b_partition)
                
                order_info = f"|A| = {q}^{a_order_m}, |B| = {q}^{b_order_m}"
                
                final_output = f"{group_desc} \t (where {order_info})"

                if final_output not in seen_structures:
                    print(final_output)
                    seen_structures.add(final_output)

if __name__ == '__main__':
    # Example usage:
    q_prime = 3
    m_natural = 3
    list_filled_groups_abelian_sylow(q_prime, m_natural)
