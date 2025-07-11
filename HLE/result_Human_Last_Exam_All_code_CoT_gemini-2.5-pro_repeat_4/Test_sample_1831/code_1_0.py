import math

def solve_endomorphism_classes():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.

    This is equivalent to counting the number of non-isomorphic functional graphs
    on 4 vertices.
    """

    # Step 1: Pre-computed numbers of non-isomorphic rooted trees for k=1 to 4.
    # t[k] is the number of rooted trees on k vertices. (OEIS A000081)
    num_rooted_trees = {1: 1, 2: 1, 3: 2, 4: 4}

    # Step 2: Calculate C_k, the number of connected functional graphs on k vertices.
    # A connected functional graph has exactly one cycle. We sum over possible cycle lengths.
    C = {}

    # C_1: A single vertex with a self-loop (cycle of length 1).
    C[1] = num_rooted_trees[1]  # A rooted tree of size 1.
    
    # C_2:
    # - Cycle length 1: A rooted tree of size 2.
    # - Cycle length 2: A 2-cycle.
    C[2] = num_rooted_trees[2] + 1
    
    # C_3:
    # - Cycle length 1: Rooted trees of size 3.
    # - Cycle length 2: A 2-cycle with a rooted tree of size 1 attached.
    # - Cycle length 3: A 3-cycle.
    C[3] = num_rooted_trees[3] + 1 + 1

    # C_4:
    # - Cycle length 1: Rooted trees of size 4.
    # - Cycle length 2: A 2-cycle with trees of total size 2 attached.
    #   The 2 vertices can be from one tree of size 2 (1 way) or two trees
    #   of size 1 (2 ways: both on one cycle vertex, or one on each). Total = 3 ways.
    # - Cycle length 3: A 3-cycle with a tree of size 1 attached.
    # - Cycle length 4: A 4-cycle.
    c4_cycle2_ways = 1 + 2 # (one t_2 tree) + (two t_1 trees)
    C[4] = num_rooted_trees[4] + c4_cycle2_ways + 1 + 1

    # Step 3: Calculate the total number of graphs by summing over partitions of 4.
    
    # Partition [4]: One component of size 4.
    count_p4 = C[4]
    
    # Partition [3, 1]: One component of size 3 and one of size 1.
    count_p31 = C[3] * C[1]
    
    # Partition [2, 2]: Two components of size 2. We can choose 2 from C[2] types
    # with replacement. This is C[2]*(C[2]+1)/2.
    count_p22 = C[2] * (C[2] + 1) // 2
    
    # Partition [2, 1, 1]: One component of size 2 and two of size 1.
    # Since there's only one type of component of size 1, this is just C[2].
    count_p211 = C[2]
    
    # Partition [1, 1, 1, 1]: Four components of size 1. Only one way.
    count_p1111 = 1
    
    partitions = {
        "[4]": count_p4,
        "[3, 1]": count_p31,
        "[2, 2]": count_p22,
        "[2, 1, 1]": count_p211,
        "[1, 1, 1, 1]": count_p1111,
    }
    
    print("The number of structures for each partition of 4 are:")
    for part, count in partitions.items():
        print(f"Partition {part}: {count}")
    
    total = sum(partitions.values())
    
    # Step 4: Print the final sum calculation.
    equation_parts = [str(v) for v in partitions.values()]
    equation = " + ".join(equation_parts)
    
    print("\nThe total number of elements of E represented by F is the sum of these values:")
    print(f"{equation} = {total}")

solve_endomorphism_classes()
<<<19>>>