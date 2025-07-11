import math

def solve_del_pezzo_problem():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.

    This problem is equivalent to counting quintic etale algebras over Q whose
    discriminant is a power of 2. Such an algebra is a direct product of
    number fields, each unramified outside the prime 2.

    The calculation is performed by considering the partitions of 5, using the
    known counts (N_n) of number fields of degree n unramified outside 2.
    """

    # N_n = number of number fields of degree n unramified outside the prime 2.
    # These are established results from number theory databases.
    N = {
        1: 1,  # The field Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 7,
        5: 0
    }

    # --- Count algebras based on partitions of 5 ---

    # Partition 5: A single quintic field K. Number of choices is N_5.
    count_5 = N[5]

    # Partition 4+1: An algebra K x Q, where K is a quartic field.
    # Number of choices for K is N_4.
    count_4_1 = N[4]

    # Partition 3+2: An algebra K x L, where K is cubic and L is quadratic.
    # Number of choices is N_3 * N_2.
    count_3_2 = N[3] * N[2]

    # Partition 3+1+1: An algebra K x Q x Q, where K is cubic.
    # Number of choices for K is N_3.
    count_3_1_1 = N[3]

    # Partition 2+2+1: An algebra K x L x Q, where K and L are distinct
    # quadratic fields. We choose 2 distinct fields from the N_2 available.
    count_2_2_1 = math.comb(N[2], 2) if N[2] >= 2 else 0

    # Partition 2+1+1+1: An algebra K x Q x Q x Q, where K is a quadratic field.
    # Number of choices for K is N_2.
    count_2_1_1_1 = N[2]

    # Partition 1+1+1+1+1: The trivial algebra Q x Q x Q x Q x Q.
    # There is only one such algebra.
    count_1_1_1_1_1 = 1

    # The total number is the sum of these counts.
    total = count_5 + count_4_1 + count_3_2 + count_3_1_1 + count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1

    print("The total number of isomorphism classes is found by summing the possibilities for each partition of 5:")
    print(f"  - Quintic field (partition 5): {count_5}")
    print(f"  - Quartic x Q (partition 4+1): {count_4_1}")
    print(f"  - Cubic x Quadratic (partition 3+2): {count_3_2}")
    print(f"  - Cubic x Q x Q (partition 3+1+1): {count_3_1_1}")
    print(f"  - Quadratic x Quadratic x Q (partition 2+2+1): {count_2_2_1}")
    print(f"  - Quadratic x Q x Q x Q (partition 2+1+1+1): {count_2_1_1_1}")
    print(f"  - Q x Q x Q x Q x Q (partition 1+1+1+1+1): {count_1_1_1_1_1}")
    print("\nThe final calculation is:")
    print(f"{count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total}")

solve_del_pezzo_problem()
<<<14>>>