import math

def solve_del_pezzo_count():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over the rational numbers with good reduction everywhere except possibly at the prime 2.

    This is equivalent to counting the number of quintic etale Q-algebras unramified
    outside of the prime 2. This count is based on the number of number fields
    of a given degree whose discriminant is a power of 2.
    """

    # N_k is the number of number fields of degree k over Q ramified only at the prime 2.
    # These values are taken from standard number theory databases (e.g., LMFDB).
    N1 = 1  # The rational numbers, Q
    N2 = 3  # The quadratic fields Q(i), Q(sqrt(2)), Q(i*sqrt(2))
    N3 = 0  # There are no such cubic fields.
    N4 = 7  # There are 7 such quartic fields.
    N5 = 0  # There are no such quintic fields.

    # We count the number of etale algebras based on the partitions of 5.
    
    # Partition 5: A single quintic field.
    # Number of choices: N5
    count_5 = N5

    # Partition 4+1: A quartic field times Q.
    # Number of choices for the quartic field: N4
    count_4_1 = N4

    # Partition 3+2: A cubic field times a quadratic field.
    # Number of choices: N3 * N2
    count_3_2 = N3 * N2
    
    # Partition 3+1+1: A cubic field times Q times Q.
    # Number of choices for the cubic field: N3
    count_3_1_1 = N3

    # Partition 2+2+1: Two distinct quadratic fields times Q.
    # This is the number of ways to choose 2 fields from the N2 available.
    count_2_2_1 = math.comb(N2, 2)

    # Partition 2+1+1+1: One quadratic field times Q times Q times Q.
    # Number of choices for the quadratic field: N2
    count_2_1_1_1 = N2

    # Partition 1+1+1+1+1: Q x Q x Q x Q x Q.
    # This corresponds to the trivial algebra. There is only one.
    count_1_1_1_1_1 = 1
    
    total_count = count_5 + count_4_1 + count_3_2 + count_3_1_1 + count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1

    print("The total number of isomorphism classes is the sum of counts for each partition of 5.")
    print("The counts for each partition are:")
    print(f"Partition 5: {count_5}")
    print(f"Partition 4+1: {count_4_1}")
    print(f"Partition 3+2: {count_3_2}")
    print(f"Partition 3+1+1: {count_3_1_1}")
    print(f"Partition 2+2+1: C({N2}, 2) = {count_2_2_1}")
    print(f"Partition 2+1+1+1: {count_2_1_1_1}")
    print(f"Partition 1+1+1+1+1: {count_1_1_1_1_1}")
    print("\nFinal Calculation:")
    print(f"{count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total_count}")

solve_del_pezzo_count()