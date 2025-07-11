import math

def find_partitions(n):
    """
    Generator function to find all integer partitions of n.
    A partition of n is a way of writing n as a sum of positive integers.
    For example, partitions of 4 are [4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1].
    
    Args:
        n (int): The integer to partition.

    Yields:
        list: A list of integers representing a partition.
    """
    def generate(target, max_val):
        if target == 0:
            yield []
            return
        
        # Recursively find partitions by subtracting numbers from the target
        # from max_val down to 1 to avoid duplicate partitions.
        for i in range(min(target, max_val), 0, -1):
            for p in generate(target - i, i):
                yield [i] + p

    # Start the generation process
    return generate(n, n)

def find_nonabelian_filled_groups(q, m):
    """
    Lists the nonabelian filled groups of order 2*q^m for a given odd prime q
    and natural number m.
    
    Args:
        q (int): An odd prime number.
        m (int): A natural number.
    """
    # --- Introduction and Explanation ---
    print(f"Finding the nonabelian filled groups of order 2 * q^m = 2 * {q}^{m} = {2 * (q**m)}.")
    print("\nBased on a theorem from group theory, a finite group G is 'filled' if and only if")
    print("it is a generalized dihedral group Dih(H) over an abelian subgroup H of odd order.")
    print("This means G = H ⋊ C₂, where the non-identity element of C₂ inverts every element of H.")
    
    print(f"\nFor a group of order 2*{q}^{m}, the subgroup H must have order {q}^{m}.")
    print(f"Since q={q} is an odd prime, the order of H is odd. This fits the theorem.")
    print("These groups are nonabelian as long as m > 0.")
    
    print("\nThe different non-isomorphic groups of this type are determined by the structure of H.")
    print(f"The possible structures for an abelian group H of order {q}^{m} correspond to the integer partitions of m={m}.")
    print("-" * 50)

    # --- Enumeration ---
    partitions = list(find_partitions(m))
    print(f"Found {len(partitions)} partitions for m = {m}, which means there are {len(partitions)} such nonabelian filled groups.\n")

    for i, p in enumerate(partitions):
        # Build the string description for the abelian group H
        h_factors = []
        for part_val in p:
            # For a partition m = m1 + m2 + ..., H = C_(q^m1) x C_(q^m2) x ...
            # The final equation requires each number, so we calculate q^part_val
            cyclic_order = q**part_val
            h_factors.append(f"C_{cyclic_order}")
        
        h_description = " x ".join(h_factors)
        
        # Output the result for this partition
        print(f"Group {i+1}:")
        print(f"  - Corresponds to partition of m={m}: {p}")
        print(f"  - Abelian subgroup H has structure: H = {h_description}")
        print(f"  - The resulting filled group G is Dih(H) = ({h_description}) ⋊ C₂")
        print()

if __name__ == '__main__':
    # --- Parameters ---
    # q must be an odd prime, m must be a natural number >= 1
    # Example: q=3, m=4. Order = 2 * 3^4 = 162
    example_q = 3
    example_m = 4

    # --- Validation (simple) ---
    if example_q % 2 == 0 or example_q < 3:
        print("Error: q must be an odd prime.")
    elif example_m < 1:
        print("Error: m must be a natural number (>= 1).")
    else:
        find_nonabelian_filled_groups(example_q, example_m)
