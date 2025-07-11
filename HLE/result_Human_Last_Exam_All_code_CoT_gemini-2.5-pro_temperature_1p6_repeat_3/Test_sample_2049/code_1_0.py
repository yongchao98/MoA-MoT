import math

def main():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over the rational numbers with good reduction everywhere except possibly at the prime 2.
    
    This is equivalent to counting quintic etale algebras over Q unramified outside 2.
    """
    
    # n[k] is the number of number fields of degree k over Q unramified outside the prime 2.
    # These values are based on standard databases and theorems in algebraic number theory.
    n = {
        1: 1,  # The field Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 6,
        5: 0,
    }

    print("The total number is found by summing the number of possible algebras for each partition of 5.")
    print("The counts for each partition are:")
    
    # Partition 5: A single quintic field.
    # The algebra is a single number field of degree 5.
    # The number of choices is n_5.
    count_5 = n[5]
    print(f"- Partition (5): {count_5}")

    # Partition 4+1: A quartic field times Q.
    # The algebra is K_4 x Q, where K_4 is a field of degree 4.
    # The number of choices is n_4 * n_1 = n_4.
    count_4_1 = n[4] * n[1]
    print(f"- Partition (4,1): {count_4_1}")

    # Partition 3+2: A cubic field times a quadratic field.
    # The algebra is K_3 x K_2.
    # The number of choices is n_3 * n_2.
    count_3_2 = n[3] * n[2]
    print(f"- Partition (3,2): {count_3_2}")
    
    # Partition 3+1+1: A cubic field times Q x Q.
    # The algebra is K_3 x Q x Q.
    # The number of choices is n_3.
    count_3_1_1 = n[3]
    print(f"- Partition (3,1,1): {count_3_1_1}")

    # Partition 2+2+1: Two quadratic fields times Q.
    # The algebra is K_2a x K_2b x Q. We are choosing a multiset of size 2
    # from the n_2 available quadratic fields.
    # The formula for multisets of size k from a set of size n is C(n+k-1, k).
    k = 2
    count_2_2_1 = math.comb(n[2] + k - 1, k)
    print(f"- Partition (2,2,1): {count_2_2_1}")

    # Partition 2+1+1+1: A quadratic field times Q x Q x Q.
    # The algebra is K_2 x Q x Q x Q.
    # The number of choices is n_2.
    count_2_1_1_1 = n[2]
    print(f"- Partition (2,1,1,1): {count_2_1_1_1}")
    
    # Partition 1+1+1+1+1: Q x Q x Q x Q x Q.
    # The "split" case. There is only one such algebra.
    count_1_1_1_1_1 = 1
    print(f"- Partition (1,1,1,1,1): {count_1_1_1_1_1}")
    
    # Summing the counts from all partitions.
    total = (count_5 + count_4_1 + count_3_2 + count_3_1_1 + 
             count_2_2_1 + count_2_1_1_1 + count_1_1_1_1_1)

    print("\nThe total number of isomorphism classes is the sum of these counts.")
    
    # Final equation as requested.
    print(f"Total = {count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1}")
    
    print(f"\nThe final answer is: {total}")

if __name__ == "__main__":
    main()
