import collections

def get_partitions(m):
    """
    Generates all unique partitions of an integer m.
    A partition is a way of writing m as a sum of positive integers.
    For example, partitions of 4 are: [4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1].
    This function uses a known concise generator algorithm and returns a list
    of unique partitions, sorted canonically.
    """
    # Using a set to store partitions ensures uniqueness after sorting.
    partitions_set = set()

    # The recursive generator for partitions.
    # From a StackOverflow post, attributed to John M. Zelle.
    def generate(n):
        # Base case: one way to partition 0 is with an empty list.
        if n == 0:
            yield []
            return
        
        # Inductive step: partitions of n are created from partitions of n-1.
        for p in generate(n - 1):
            # Option 1: Add a 1 to the partition.
            yield [1] + p
            # Option 2: Increment the largest part of the partition.
            # This avoids duplicates by only incrementing when the parts are in
            # non-decreasing order.
            if p and (len(p) < 2 or p[1] > p[0]):
                yield [p[0] + 1] + p[1:]

    for p in generate(m):
        # Add the canonical form (sorted descending) to the set.
        partitions_set.add(tuple(sorted(p, reverse=True)))

    # Return a sorted list of partitions for deterministic output.
    return sorted(list(partitions_set), reverse=True)


def describe_filled_groups(m, q_char='q'):
    """
    Based on the classification of finite filled groups, this function describes
    the nonabelian filled groups of order 2*q^m.
    
    Args:
        m (int): A natural number exponent.
        q_char (str): The character to use for the odd prime q.
    """

    print(
f"""Based on group theory, the nonabelian filled groups of order 2*{q_char}^{m}
(for an odd prime {q_char} and natural number m={m}) are the generalized dihedral groups.
These groups have the structure G = A ⋊ C₂, where:
- A is an abelian group of order {q_char}^{m}.
- C₂ is the cyclic group of order 2.
- The semidirect product action is inversion (a ↦ a⁻¹ for all a ∈ A).

The structure of the abelian group A is determined by the partitions of m.
The following are the possible structures for A and the corresponding group G:
"""
    )
    
    partitions = get_partitions(m)
    
    print(f"Found {len(partitions)} possible group structures for m = {m}.\n")
    
    for i, p in enumerate(partitions):
        # Format the partition sum equation, e.g., "4 = 2 + 2"
        partition_eq = f"{m} = {' + '.join(map(str, p))}"

        # Format the abelian group structure string, e.g., "C_(q^2) x C_(q^2)"
        a_parts = []
        for power in p:
            if power == 1:
                a_parts.append(f"C_{q_char}")
            else:
                a_parts.append(f"C_({q_char}^{power})")
        a_structure = " x ".join(a_parts)
        
        print(f"Structure {i+1}:")
        print(f"  Partition: {partition_eq}")
        print(f"  Abelian component A = {a_structure}")
        print(f"  The filled group is G = ({a_structure}) ⋊ C₂\n")


if __name__ == '__main__':
    # You can change the value of m to explore different orders.
    # For example, try m = 2, 3, 4, or 5.
    m_exponent = 4
    
    # Check if m is a valid natural number.
    if not isinstance(m_exponent, int) or m_exponent < 1:
        print("Error: m must be a natural number (positive integer).")
    else:
        describe_filled_groups(m_exponent)
