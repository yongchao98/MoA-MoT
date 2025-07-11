import math
from collections import Counter
import itertools

def main():
    """
    Calculates and lists all non-isomorphic finite Weyl groups of rank 4.
    """
    # Step 1: Define the non-isomorphic irreducible Weyl groups by rank.
    # Note: B_n and C_n are isomorphic. D_n is irreducible for n>=4.
    # A_3 is isomorphic to D_3.
    irred_by_rank = {
        1: ["A1"],
        2: ["A2", "G2"],
        3: ["A3", "B3"],
        4: ["A4", "B4", "D4", "F4"]
    }

    # The partitions of the integer 4
    partitions = [
        [4],
        [3, 1],
        [2, 2],
        [2, 1, 1],
        [1, 1, 1, 1]
    ]

    print("The non-isomorphic finite Weyl groups of rank 4 are determined by the partitions of 4.\n")
    
    total_count = 0
    counts_per_partition = []

    # Iterate through each partition of 4
    for p in partitions:
        partition_str = " + ".join(map(str, p))
        print(f"--- For Partition {partition_str} ---")

        # Count occurrences of each part, e.g., [2, 1, 1] -> {2: 1, 1: 2}
        part_counts = Counter(p)
        
        # Get the list of possible groups for each unique part
        options_per_part = {
            part: irred_by_rank.get(part, []) for part in part_counts
        }

        # Generate combinations for this partition
        # For each unique part, we generate combinations with replacement
        # E.g., for part 2 (count 2), from [A2, G2], we can have (A2,A2), (A2,G2), (G2,G2)
        group_choices = []
        for part, count in part_counts.items():
            choices = list(itertools.combinations_with_replacement(options_per_part[part], count))
            group_choices.append(choices)

        # Get the cartesian product of the choices for each part
        # E.g., for partition [3, 1], we have choices for 3: ((A3,), (B3,))
        # and choices for 1: ((A1),). The product is ((A3,),(A1,)), ((B3,),(A1,))
        count_for_this_partition = 0
        current_groups = []
        for combination in itertools.product(*group_choices):
            # Flatten the list of tuples into a single list of group names
            full_group_list = [item for sublist in combination for item in sublist]
            group_name = " x ".join(sorted(full_group_list))
            current_groups.append(group_name)
            count_for_this_partition += 1
        
        for name in sorted(current_groups):
            print(name)

        print(f"Count for this partition: {count_for_this_partition}\n")
        counts_per_partition.append(str(count_for_this_partition))
        total_count += count_for_this_partition

    # Final summary
    print("------------------------------")
    print("The total number is the sum of counts from each partition.")
    final_equation = " + ".join(counts_per_partition)
    print(f"Total = {final_equation} = {total_count}")


if __name__ == "__main__":
    main()
