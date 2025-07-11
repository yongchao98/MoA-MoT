import math

def main():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction outside the prime 2.

    This is equivalent to counting the number of quartic etale algebras over Q
    unramified outside 2.
    """

    # Number of number fields unramified outside the prime 2 for each degree n <= 4.
    # This data comes from number theory databases (e.g., LMFDB).
    fields_count = {
        1: 1,  # The field Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 7,
    }

    print("This program calculates the number of del Pezzo fibrations of degree 5 over Spec(Z[1/2]).")
    print("This is equivalent to counting quartic etale algebras over Q unramified outside prime 2.")
    print("\nFirst, we use the known counts of number fields unramified outside 2:")
    for degree, count in fields_count.items():
        print(f"   - Degree {degree}: {count} field(s)")

    print("\nNext, we count the algebras for each partition of the total degree 4:")

    # Partition 4: A single quartic field
    # This corresponds to an irreducible quartic etale algebra.
    count_4 = fields_count[4]
    print(f"Partition 4: {count_4} (from quartic fields)")

    # Partition 3+1: A cubic field x a linear field
    count_3_1 = fields_count[3] * fields_count[1]
    print(f"Partition 3+1: {fields_count[3]} * {fields_count[1]} = {count_3_1}")

    # Partition 2+2: A product of two quadratic fields
    # This is choosing 2 fields from the set of 3 quadratic fields, with replacement.
    # The formula for multisets of size k from a set of size n is C(n+k-1, k).
    n = fields_count[2]
    k = 2
    # C(3+2-1, 2) = C(4, 2) = 6
    count_2_2 = math.comb(n + k - 1, k)
    print(f"Partition 2+2: {count_2_2} (product of two quadratic fields, chosen from {n} with replacement)")

    # Partition 2+1+1: A quadratic field x linear x linear
    # The linear field is just Q, so we only need to choose the quadratic field.
    count_2_1_1 = fields_count[2]
    print(f"Partition 2+1+1: {count_2_1_1} (one quadratic field and two copies of Q)")

    # Partition 1+1+1+1: Four copies of the linear field Q
    count_1_1_1_1 = 1
    print(f"Partition 1+1+1+1: {count_1_1_1_1} (four copies of Q)")

    # Summing up all the possibilities
    total_count = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1

    print("\nThe total number is the sum of the counts for each partition.")
    print(f"Total = {count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_count}")

if __name__ == "__main__":
    main()
