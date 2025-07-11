import math

def main():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction outside the prime 2.

    This is equivalent to counting the number of quintic étale algebras over Q
    unramified outside 2.
    """

    # n_d is the number of number fields of degree d unramified outside the prime 2.
    # These values are obtained from number theory databases like LMFDB.
    n = {
        1: 1,  # Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 7,
        5: 6
    }

    # We count the number of quintic etale algebras by considering the partitions of 5.
    # An etale algebra L is a product of number fields K_i where sum([K_i:Q]) = 5.

    # Partition 5: A single quintic field.
    # L = K_5
    count_5 = n[5]

    # Partition 4+1: A quartic field and a degree 1 field (Q).
    # L = K_4 x Q
    count_4_1 = n[4]

    # Partition 3+2: A cubic field and a quadratic field.
    # L = K_3 x K_2
    count_3_2 = n[3] * n[2]

    # Partition 3+1+1: A cubic field and Q x Q.
    # L = K_3 x Q x Q
    count_3_1_1 = n[3]

    # Partition 2+2+1: Two quadratic fields and Q.
    # L = K_{2,a} x K_{2,b} x Q
    # This is equivalent to choosing 2 fields from a set of n[2]=3 with replacement.
    # The number of multisets of size k from a set of size n is C(n+k-1, k).
    count_2_2_1 = math.comb(n[2] + 2 - 1, 2)

    # Partition 2+1+1+1: One quadratic field and Q x Q x Q.
    # L = K_2 x Q x Q x Q
    count_2_1_1_1 = n[2]

    # Partition 1+1+1+1+1: Five copies of Q.
    # L = Q x Q x Q x Q x Q
    count_1_1_1_1_1 = 1

    counts = {
        "5": count_5,
        "4+1": count_4_1,
        "3+2": count_3_2,
        "3+1+1": count_3_1_1,
        "2+2+1": count_2_2_1,
        "2+1+1+1": count_2_1_1_1,
        "1+1+1+1+1": count_1_1_1_1_1
    }
    
    print("The number of isomorphism classes is the sum of counts for each partition of 5:")
    
    total = 0
    equation_parts = []
    for partition, count in counts.items():
        print(f"Partition {partition}: {count}")
        total += count
        equation_parts.append(str(count))
        
    equation = " + ".join(equation_parts)
    print(f"\nFinal calculation: {equation} = {total}")

if __name__ == "__main__":
    main()
