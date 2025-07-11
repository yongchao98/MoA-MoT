import math

def calculate_del_pezzo_classes():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.

    This is equivalent to counting the number of etale Q-algebras of degree 5
    unramified outside the prime 2.
    """

    # Step 1: Define the number of number fields of degree k unramified outside the prime 2.
    # This data is from the L-Functions and Modular Forms Database (LMFDB).
    N = {
        1: 1,  # Q
        2: 3,  # Q(sqrt(-1)), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 6,
        5: 0,
    }

    print("This script calculates the number of del Pezzo fibrations of degree 5 over Spec(Z[1/2]).")
    print("This is equivalent to counting etale Q-algebras of degree 5 unramified outside {2}.")
    print("\nFirst, we list the number of number fields of a given degree unramified outside {2}:")
    for degree, count in N.items():
        print(f"  - Degree {degree}: {count} fields")

    print("\nNext, we sum the number of possible algebras for each partition of 5:")

    # We will store the results for each partition to display in the final sum.
    partition_counts = {}

    # Partition 5
    # Algebra is a single field of degree 5.
    count_5 = N[5]
    partition_counts["5"] = count_5

    # Partition 4+1
    # Algebra is a product of a degree 4 field and a degree 1 field.
    count_4_1 = N[4] * N[1]
    partition_counts["4+1"] = count_4_1

    # Partition 3+2
    # Algebra is a product of a degree 3 field and a degree 2 field.
    count_3_2 = N[3] * N[2]
    partition_counts["3+2"] = count_3_2
    
    # Partition 3+1+1
    # Algebra is a product of a degree 3 field and two degree 1 fields.
    # Choose 1 from N[3] and 2 from N[1] with replacement.
    count_3_1_1 = N[3] * math.comb(N[1] + 2 - 1, 2)
    partition_counts["3+1+1"] = count_3_1_1

    # Partition 2+2+1
    # Algebra is a product of two degree 2 fields and one degree 1 field.
    # Choose 2 from N[2] with replacement and 1 from N[1].
    count_2_2_1 = math.comb(N[2] + 2 - 1, 2) * N[1]
    partition_counts["2+2+1"] = count_2_2_1

    # Partition 2+1+1+1
    # Algebra is a product of one degree 2 field and three degree 1 fields.
    # Choose 1 from N[2] and 3 from N[1] with replacement.
    count_2_1_1_1 = N[2] * math.comb(N[1] + 3 - 1, 3)
    partition_counts["2+1+1+1"] = count_2_1_1_1

    # Partition 1+1+1+1+1
    # Algebra is a product of five degree 1 fields.
    # Choose 5 from N[1] with replacement.
    count_1_1_1_1_1 = math.comb(N[1] + 5 - 1, 5)
    partition_counts["1+1+1+1+1"] = count_1_1_1_1_1

    total_count = sum(partition_counts.values())

    # Print the contribution of each partition
    for part, count in partition_counts.items():
        print(f"  - Partition {part.ljust(11)}: {count} classes")

    # Final result construction
    final_sum_str = " + ".join(str(c) for c in partition_counts.values())
    print("\nFinal Calculation:")
    print(f"Total number = {final_sum_str} = {total_count}")

calculate_del_pezzo_classes()