import math

def calculate_and_print_classes():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.
    This is equivalent to counting non-isomorphic functional graphs on 4 vertices.
    The method is to sum the number of possibilities for each integer partition of 4.
    """
    print("The total number of classes is found by summing the number of possibilities for each partition of the set of 4 elements.\n")

    # Case 1: Partition [4] (a single component of 4 vertices)
    # The structure is determined by the length of the cycle within the component.
    # - Cycle of length 4: 1 type (a simple 4-cycle).
    # - Cycle of length 3: 1 type (a 3-cycle with one vertex attached).
    # - Cycle of length 2: 3 types (based on how the two non-cycle vertices attach).
    # - Cycle of length 1: 4 types (the number of non-isomorphic rooted trees on 4 vertices).
    count_part_4 = 1 + 1 + 3 + 4
    print(f"For partition [4], number of classes = 1 (4-cycle) + 1 (3-cycle) + 3 (2-cycle) + 4 (1-cycle) = {count_part_4}")

    # Case 2: Partition [3, 1] (one component of size 3, one of size 1)
    # Number of structures for a size 3 component:
    # - Cycle of length 3: 1 type
    # - Cycle of length 2: 1 type
    # - Cycle of length 1 (rooted tree on 3 vertices): 2 types
    # Total for size 3 component = 1 + 1 + 2 = 4
    # Number of structures for a size 1 component:
    # - Cycle of length 1 (a fixed point): 1 type
    count_part_3_1 = (1 + 1 + 2) * 1
    print(f"For partition [3, 1], number of classes = (1 + 1 + 2) * 1 = {count_part_3_1}")

    # Case 3: Partition [2, 2] (two components of size 2)
    # Number of structures for a size 2 component:
    # - A 2-cycle: 1 type
    # - A rooted tree on 2 vertices: 1 type
    # Total types for a size 2 component = 2.
    # We need to choose 2 components from these 2 types, with replacement.
    # The combinations are {2-cycle, 2-cycle}, {2-cycle, tree}, {tree, tree}. This gives 3 ways.
    count_part_2_2 = 3
    print(f"For partition [2, 2], number of classes = {count_part_2_2}")

    # Case 4: Partition [2, 1, 1] (one component of size 2, two of size 1)
    # Number of types for size 2 component = 2.
    # Number of types for size 1 component = 1.
    # The two size 1 components are of the same type. So, we just multiply the number of types.
    count_part_2_1_1 = 2 * 1
    print(f"For partition [2, 1, 1], number of classes = 2 * 1 = {count_part_2_1_1}")

    # Case 5: Partition [1, 1, 1, 1] (four components of size 1)
    # Each component is a fixed point. There is only one way to have four fixed points.
    count_part_1_1_1_1 = 1
    print(f"For partition [1, 1, 1, 1], number of classes = {count_part_1_1_1_1}")

    # Final Summation
    total_count = count_part_4 + count_part_3_1 + count_part_2_2 + count_part_2_1_1 + count_part_1_1_1_1
    print("\nTotal number of equivalence classes is the sum of the counts for each partition:")
    print(f"{count_part_4} + {count_part_3_1} + {count_part_2_2} + {count_part_2_1_1} + {count_part_1_1_1_1} = {total_count}")

calculate_and_print_classes()