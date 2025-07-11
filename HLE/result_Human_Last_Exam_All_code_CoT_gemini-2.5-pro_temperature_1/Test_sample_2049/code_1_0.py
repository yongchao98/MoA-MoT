import math

def solve_del_pezzo_problem():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.

    This is equivalent to counting the number of quintic etale algebras over Q
    unramified outside the prime 2.
    """

    # Step 1: Number of number fields of degree n unramified outside prime 2.
    # These are known results from number theory databases (e.g., LMFDB).
    A = {
        1: 1,  # The field Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,  # No such cubic fields exist
        4: 7,  # 7 such quartic fields exist
        5: 0   # No such quintic fields exist
    }

    # Step 2: Count the number of quintic etale algebras for each partition of 5.
    # An etale algebra is a product of distinct field extensions.
    # Let's calculate the number of algebras for each partition of 5.

    # Partition 5: A single quintic field.
    # L5
    p_5 = A[5]

    # Partition 4+1: A quartic field times Q.
    # L4 x Q
    p_4_1 = A[4] * A[1]

    # Partition 3+2: A cubic field times a quadratic field.
    # L3 x L2
    p_3_2 = A[3] * A[2]
    
    # Partition 3+1+1: A cubic field times Q x Q.
    # L3 x Q x Q
    p_3_1_1 = A[3]

    # Partition 2+2+1: Two distinct quadratic fields times Q.
    # L2_a x L2_b x Q
    # We choose 2 distinct fields from the A[2]=3 available quadratic fields.
    p_2_2_1 = math.comb(A[2], 2)

    # Partition 2+1+1+1: One quadratic field times Q x Q x Q.
    # L2 x Q x Q x Q
    p_2_1_1_1 = A[2]

    # Partition 1+1+1+1+1: Q x Q x Q x Q x Q.
    # There is only one such algebra.
    p_1_1_1_1_1 = 1

    # Step 3: Sum the counts for all partitions.
    total = p_5 + p_4_1 + p_3_2 + p_3_1_1 + p_2_2_1 + p_2_1_1_1 + p_1_1_1_1_1
    
    # Print the equation with each term.
    print("The total number is the sum of possibilities from each partition of 5:")
    print(f"Partition (5): {p_5}")
    print(f"Partition (4,1): {p_4_1}")
    print(f"Partition (3,2): {p_3_2}")
    print(f"Partition (3,1,1): {p_3_1_1}")
    print(f"Partition (2,2,1): {p_2_2_1}")
    print(f"Partition (2,1,1,1): {p_2_1_1_1}")
    print(f"Partition (1,1,1,1,1): {p_1_1_1_1_1}")
    print("\nFinal equation:")
    print(f"{p_5} + {p_4_1} + {p_3_2} + {p_3_1_1} + {p_2_2_1} + {p_2_1_1_1} + {p_1_1_1_1_1} = {total}")

    print(f"\nThe total number is: {total}")

solve_del_pezzo_problem()
<<<14>>>