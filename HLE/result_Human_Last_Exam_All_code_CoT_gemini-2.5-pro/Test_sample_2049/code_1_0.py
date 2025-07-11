import math

def solve():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.

    This is equivalent to counting the number of quintic etale algebras over Q
    unramified outside {2, infinity}.
    """

    # f[n] stores the number of number fields of degree n over Q
    # unramified outside {2, infinity}. These values are from number theory databases (LMFDB).
    f = {
        1: 1,  # Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 9,
        5: 6
    }

    # We calculate the number of etale algebras for each partition of 5.

    # Partition 5: A single number field of degree 5.
    # Number of choices is f[5].
    count_5 = f[5]

    # Partition 4+1: A product of a degree 4 field and a degree 1 field.
    # Number of choices is f[4] * f[1].
    count_4_1 = f[4] * f[1]

    # Partition 3+2: A product of a degree 3 field and a degree 2 field.
    # Number of choices is f[3] * f[2].
    count_3_2 = f[3] * f[2]

    # Partition 3+1+1: A product of a degree 3 field and two degree 1 fields.
    # Since f[3] is 0, this is 0.
    # The number of ways to choose 2 fields of degree 1 with replacement is math.comb(f[1] + 2 - 1, 2)
    count_3_1_1 = f[3] * math.comb(f[1] + 2 - 1, 2)

    # Partition 2+2+1: A product of two degree 2 fields and one degree 1 field.
    # We choose 2 fields of degree 2 with replacement from f[2] options.
    count_2_2_1 = math.comb(f[2] + 2 - 1, 2) * f[1]

    # Partition 2+1+1+1: A product of one degree 2 field and three degree 1 fields.
    # We choose 1 field of degree 2 from f[2] options.
    count_2_1_1_1 = f[2] * math.comb(f[1] + 3 - 1, 3)

    # Partition 1+1+1+1+1: A product of five degree 1 fields.
    # We choose 5 fields of degree 1 with replacement from f[1] option.
    count_1_1_1_1_1 = math.comb(f[1] + 5 - 1, 5)

    total_count = (count_5 + count_4_1 + count_3_2 + count_3_1_1 +
                   count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1)

    print("The total number is the sum of counts for each partition of 5:")
    print(f"Partition 5: {count_5}")
    print(f"Partition 4+1: {count_4_1}")
    print(f"Partition 3+2: {count_3_2}")
    print(f"Partition 3+1+1: {count_3_1_1}")
    print(f"Partition 2+2+1: {count_2_2_1}")
    print(f"Partition 2+1+1+1: {count_2_1_1_1}")
    print(f"Partition 1+1+1+1+1: {count_1_1_1_1_1}")
    print("-" * 20)
    print("The final equation is:")
    print(f"{count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total_count}")
    print("-" * 20)
    print(f"The total number of isomorphism classes is {total_count}.")


if __name__ == '__main__':
    solve()