import math

def main():
    """
    Calculates and lists the number of non-isomorphic finite Weyl groups of rank 4.
    """
    # Define the irreducible Weyl groups by their rank up to rank 4.
    # We use standard representatives to handle isomorphisms, e.g., B_n for C_n, and A_3 for D_3.
    groups_by_rank = {
        1: ["A_1"],
        2: ["A_2", "B_2", "G_2"],
        3: ["A_3", "B_3"],
        4: ["A_4", "B_4", "D_4", "F_4"]
    }

    all_groups = []
    counts_per_case = []

    # Case 1: One component of rank 4
    # Partition: 4
    case1_groups = list(groups_by_rank[4])
    all_groups.extend(case1_groups)
    counts_per_case.append(len(case1_groups))

    # Case 2: Two components with ranks summing to 4
    # Partition: 3 + 1
    case2a_groups = []
    for g3 in groups_by_rank[3]:
        for g1 in groups_by_rank[1]:
            case2a_groups.append(f"{g3} x {g1}")
    
    # Partition: 2 + 2
    case2b_groups = []
    rank2_list = groups_by_rank[2]
    # To find non-isomorphic combinations (order doesn't matter), we iterate such that j >= i.
    # This is equivalent to combinations with replacement.
    for i in range(len(rank2_list)):
        for j in range(i, len(rank2_list)):
            case2b_groups.append(f"{rank2_list[i]} x {rank2_list[j]}")
            
    all_groups.extend(case2a_groups)
    all_groups.extend(case2b_groups)
    # The count for 2 components is the sum of the two sub-cases
    counts_per_case.append(len(case2a_groups))
    counts_per_case.append(len(case2b_groups))


    # Case 3: Three components with ranks summing to 4
    # Partition: 2 + 1 + 1
    case3_groups = []
    # The two rank 1 components are always A_1 x A_1
    for g2 in groups_by_rank[2]:
        case3_groups.append(f"{g2} x {groups_by_rank[1][0]} x {groups_by_rank[1][0]}")
    all_groups.extend(case3_groups)
    counts_per_case.append(len(case3_groups))

    # Case 4: Four components with ranks summing to 4
    # Partition: 1 + 1 + 1 + 1
    g1 = groups_by_rank[1][0]
    case4_groups = [f"{g1} x {g1} x {g1} x {g1}"]
    all_groups.extend(case4_groups)
    counts_per_case.append(len(case4_groups))

    print("The non-isomorphic finite Weyl groups of rank 4 are:")
    for group in sorted(all_groups):
        print(group)
    
    total_count = sum(counts_per_case)
    
    print("\nTotal number calculation:")
    # The prompt requires outputting each number in the final equation.
    # counts_per_case currently is [4, 2, 6, 3, 1]
    # This represents partitions (4), (3+1), (2+2), (2+1+1), (1+1+1+1)
    equation_str = f"{counts_per_case[0]} (rank 4) + {counts_per_case[1]} (rank 3+1) + {counts_per_case[2]} (rank 2+2) + {counts_per_case[3]} (rank 2+1+1) + {counts_per_case[4]} (rank 1+1+1+1)"
    
    # Let's print the final equation numbers clearly
    final_equation_numbers = " + ".join(map(str, counts_per_case))
    print(f"{final_equation_numbers} = {total_count}")

if __name__ == "__main__":
    main()
