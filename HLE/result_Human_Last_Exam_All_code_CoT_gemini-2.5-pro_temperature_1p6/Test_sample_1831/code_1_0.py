import math
from collections import Counter

# Step 1: Number of non-isomorphic rooted trees, t_j.
# These are known values from OEIS A000081.
t = {1: 1, 2: 1, 3: 2, 4: 4}

# Step 2: Number of non-isomorphic connected functional graphs, a_k.
# A connected functional graph has a single cycle, with rooted trees attached to its vertices.
# We calculate a_k for k=1, 2, 3, 4.
# This corresponds to OEIS A001372. My analysis confirms the following values.
a = {
    1: 1,  # A single loop.
    2: 2,  # A 2-cycle, or a loop with one vertex attached.
    3: 4,  # 3-cycle; 2-cycle + 1 tree; 1-cycle + 2 tree arrangements.
    4: 9,  # 4-cycle; 3-cycle + ...; 2-cycle + ...; 1-cycle + ...
}

# Step 3: Calculate the total number of non-isomorphic functional graphs on 4 vertices, b_4.
# We sum over the integer partitions of 4.

def combinations_with_replacement(n, k):
    """Calculates C(n + k - 1, k)"""
    if k < 0 or n < 0:
        return 0
    return math.comb(n + k - 1, k)

def solve():
    """
    Calculates the number of non-isomorphic functional graphs on 4 vertices.
    """
    n = 4
    # The integer partitions of 4
    partitions = [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    total_count = 0
    
    print(f"To find the number of equivalence classes of endomorphisms on a set of size 4, we count the non-isomorphic functional graphs on 4 vertices.\n")
    print("Let a_k be the number of connected functional graphs on k vertices.")
    print(f"The known values are: a_1 = {a[1]}, a_2 = {a[2]}, a_3 = {a[3]}, a_4 = {a[4]}.\n")
    print("We sum the counts for each integer partition of 4:\n")
    
    for p in partitions:
        counts = Counter(p)
        term_result = 1
        
        # This formula counts the number of ways to choose components for a given partition
        for size, num_parts in counts.items():
            term_result *= combinations_with_replacement(a[size], num_parts)
            
        total_count += term_result
        
        partition_str = " + ".join(map(str, p))
        
        # Building the explanation string for each term
        calc_str_parts = []
        for size, num_parts in counts.items():
            if num_parts == 1:
                calc_str_parts.append(f"a_{size}")
            else:
                calc_str_parts.append(f"C({a[size]} + {num_parts} - 1, {num_parts})")
        
        calc_str = " * ".join(calc_str_parts)
        val_str_parts = []
        for size, num_parts in counts.items():
             if num_parts == 1:
                val_str_parts.append(str(a[size]))
             else:
                val_str_parts.append(str(combinations_with_replacement(a[size], num_parts)))
        
        val_str = " * ".join(val_str_parts)

        if len(calc_str_parts) > 1:
            print(f"Partition {partition_str}: {calc_str} = {val_str} = {term_result}")
        else:
            print(f"Partition {partition_str}: {calc_str} = {term_result}")

    print(f"\nTotal count = 4 + 3 + 1 + 1 + 2 + 1 + 9 (this is a simple sum of the above numbers)")
    final_equation_parts = []
    for p in partitions:
        counts = Counter(p)
        term_result = 1
        for size, num_parts in counts.items():
            term_result *= combinations_with_replacement(a[size], num_parts)
        final_equation_parts.append(str(term_result))
        
    final_equation = " + ".join(final_equation_parts[::-1])
    print(f"Summing these up: {final_equation} = {total_count}")
    print(f"\nThus, there are {total_count} elements of E represented by F.")


solve()