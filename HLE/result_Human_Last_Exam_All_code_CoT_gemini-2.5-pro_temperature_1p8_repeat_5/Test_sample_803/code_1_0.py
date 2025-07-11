import collections

def get_partitions(n):
    """
    Generates all integer partitions of n.
    Yields each partition as a list of integers.
    """
    if n == 0:
        yield []
        return
    if n < 0:
        return
        
    for p in get_partitions(n - 1):
        yield [1] + p
        if p and (len(p) < 2 or p[1] > p[0]):
            yield [p[0] + 1] + p[1:]

def find_nonabelian_filled_groups(q, m):
    """
    Based on group theory results, finds and describes the nonabelian filled 
    groups of order 2 * q**m for a given odd prime q and natural number m.

    The theory indicates these are the generalized dihedral groups Dih(Q) 
    where Q is any abelian group of order q**m.
    """
    
    print(f"Finding nonabelian filled groups of order 2 * {q}^{m} = {2 * q**m}")
    print("=" * 60)
    print("The nonabelian filled groups of this order are the generalized dihedral groups Dih(Q),")
    print("where Q is an abelian group of order q^m and the action is by inversion.\n")
    print("The possible structures for Q correspond to the integer partitions of m:")
    print("-" * 60)

    # Use a set to store unique sorted partitions
    partitions = set()
    for p in get_partitions(m):
        partitions.add(tuple(sorted(p, reverse=True)))

    if not partitions:
        print("No partitions found for m=0.")
        return

    for i, p in enumerate(sorted(list(partitions), reverse=True)):
        group_q_parts = []
        for lambda_val in p:
            if q == 1: # Should not happen based on prompt (q is odd prime)
                group_q_parts.append(f"C_{lambda_val}")
            else:
                group_q_parts.append(f"C_{{{q}^{lambda_val}}}")

        group_q_str = " x ".join(group_q_parts)
        
        print(f"Group {i+1}:")
        print(f"  - This corresponds to the partition m = {' + '.join(map(str, p))}.")
        print(f"  - The abelian group Q is: {group_q_str}")
        print(f"  - The nonabelian filled group is G = Dih(Q) = ({group_q_str}) \u22ca C_2")
        print()

if __name__ == '__main__':
    # Example: q=3 (an odd prime), m=4 (a natural number)
    # The groups will have order 2 * 3^4 = 162.
    example_q = 3
    example_m = 4
    
    find_nonabelian_filled_groups(example_q, example_m)
