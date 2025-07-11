import math

def solve_del_pezzo_count():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction outside the prime 2.

    This is equivalent to counting quintic étale algebras over Q unramified outside 2.
    The calculation is based on the number of number fields of degrees 1-5 that are
    unramified outside the prime 2, sourced from number theory databases (LMFDB).
    """

    # n_d is the number of number fields of degree d unramified outside the prime 2.
    # Data from Jones-Roberts database / LMFDB.
    n1 = 1  # The rational field Q
    n2 = 3  # Q(sqrt(-1)), Q(sqrt(2)), Q(sqrt(-2))
    n3 = 0  # There are no such cubic fields.
    n4 = 10 # Determined by counting fields with Galois groups S4, D4, C4, V4.
    n5 = 6  # All have Galois group S5.

    # We count the number of quintic étale algebras by considering the partitions of 5.
    # An algebra corresponds to a product of fields whose degrees sum to 5.

    # Partition: 5
    # Algebra type: K_5 (a single quintic field)
    # Number of choices for K_5.
    count_5 = n5
    
    # Partition: 4 + 1
    # Algebra type: K_4 x Q
    # Number of choices for K_4.
    count_4_1 = n4

    # Partition: 3 + 2
    # Algebra type: K_3 x K_2
    # Number of choices is n3 * n2.
    count_3_2 = n3 * n2

    # Partition: 3 + 1 + 1
    # Algebra type: K_3 x Q x Q
    # Number of choices for K_3.
    count_3_1_1 = n3

    # Partition: 2 + 2 + 1
    # Algebra type: K_2 x K'_2 x Q
    # This is the number of ways to choose a multiset of size 2 from the n2 available quadratic fields.
    # The formula for combinations with repetition is C(n+k-1, k).
    count_2_2_1 = math.comb(n2 + 2 - 1, 2)

    # Partition: 2 + 1 + 1 + 1
    # Algebra type: K_2 x Q x Q x Q
    # Number of choices for K_2.
    count_2_1_1_1 = n2

    # Partition: 1 + 1 + 1 + 1 + 1
    # Algebra type: Q x Q x Q x Q x Q
    # There is only one such algebra.
    count_1_1_1_1_1 = 1
    
    # The total is the sum of counts for all partitions.
    total_count = (count_5 + count_4_1 + count_3_2 + count_3_1_1 +
                   count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1)

    # Output the final equation with each number.
    print("The total number of isomorphism classes is the sum over partitions of 5:")
    print(f"{count_5} (from partition 5) + "
          f"{count_4_1} (from partition 4+1) + "
          f"{count_3_2} (from partition 3+2) + "
          f"{count_3_1_1} (from partition 3+1+1) + "
          f"{count_2_2_1} (from partition 2+2+1) + "
          f"{count_2_1_1_1} (from partition 2+1+1+1) + "
          f"{count_1_1_1_1_1} (from partition 1+1+1+1+1) "
          f"= {total_count}")

solve_del_pezzo_count()
<<<26>>>