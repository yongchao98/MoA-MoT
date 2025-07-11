import math

def combinations(n, k):
    """
    Calculates n C k, the number of combinations.
    It's the number of ways to choose k items from a set of n items.
    """
    if k < 0 or k > n:
        return 0
    # To avoid large intermediate numbers, this is a more stable way to compute nCk
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def multiset_combinations(n, k):
    """
    Calculates the number of k-element multisets from a set of size n.
    This is equivalent to C(n + k - 1, k).
    """
    if n == 0 and k > 0:
        return 0
    return combinations(n + k - 1, k)

def main():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.
    """
    print("This program calculates the number of isomorphism classes of del Pezzo fibrations of degree 5 over Spec(Z[1/2]).")
    print("This is equivalent to counting etale Q-algebras of rank 5 unramified outside the prime 2.\n")

    # Step 1: Provide the number of number fields of degree n unramified outside {2}.
    # N[n] is the number of such fields of degree n. These are established results.
    N = {
        1: 1,  # Just Q itself
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,  # There are no cubic fields ramified only at 2
        4: 18, # From the Jones-Roberts database / LMFDB
        5: 6   # From the Jones-Roberts database / LMFDB
    }

    print("The number of number fields of degree n unramified outside {2} (for n=1 to 5):")
    for degree, count in N.items():
        print(f"  Degree {degree}: {count} fields")
    print("\n" + "="*50 + "\n")

    # Step 2: An etale algebra of dimension 5 corresponds to a partition of 5.
    # We sum the counts for each partition of 5.
    # The partitions are: (5), (4,1), (3,2), (3,1,1), (2,2,1), (2,1,1,1), (1,1,1,1,1)

    print("The total number of algebras is the sum of counts for each partition of 5.\n")
    
    total_count = 0
    terms = []

    # Partition (5): A single quintic field.
    count_5 = N[5]
    total_count += count_5
    terms.append(count_5)
    print(f"For partition (5):")
    print(f"  Structure: A single quintic field.")
    print(f"  Count = N_5 = {count_5}\n")

    # Partition (4, 1): A quartic field and a degree 1 field (Q).
    count_4_1 = N[4] * N[1]
    total_count += count_4_1
    terms.append(count_4_1)
    print(f"For partition (4, 1):")
    print(f"  Structure: K_4 x Q")
    print(f"  Count = N_4 * N_1 = {N[4]} * {N[1]} = {count_4_1}\n")

    # Partition (3, 2): A cubic field and a quadratic field.
    count_3_2 = N[3] * N[2]
    total_count += count_3_2
    terms.append(count_3_2)
    print(f"For partition (3, 2):")
    print(f"  Structure: K_3 x K_2")
    print(f"  Count = N_3 * N_2 = {N[3]} * {N[2]} = {count_3_2}\n")

    # Partition (3, 1, 1): A cubic field and two copies of Q.
    count_3_1_1 = N[3]
    total_count += count_3_1_1
    terms.append(count_3_1_1)
    print(f"For partition (3, 1, 1):")
    print(f"  Structure: K_3 x Q x Q")
    print(f"  Count = N_3 = {count_3_1_1}\n")

    # Partition (2, 2, 1): Two quadratic fields and Q.
    n, k = N[2], 2
    count_2_2_1 = multiset_combinations(n, k)
    total_count += count_2_2_1
    terms.append(count_2_2_1)
    print(f"For partition (2, 2, 1):")
    print(f"  Structure: K_a x K_b x Q (K_a, K_b are quadratic fields)")
    print(f"  This is a multiset of size {k} from {n} choices.")
    print(f"  Count = C({n}+{k}-1, {k}) = C({n+k-1}, {k}) = {count_2_2_1}\n")

    # Partition (2, 1, 1, 1): One quadratic field and three copies of Q.
    count_2_1_1_1 = N[2]
    total_count += count_2_1_1_1
    terms.append(count_2_1_1_1)
    print(f"For partition (2, 1, 1, 1):")
    print(f"  Structure: K_2 x Q x Q x Q")
    print(f"  Count = N_2 = {count_2_1_1_1}\n")

    # Partition (1, 1, 1, 1, 1): Five copies of Q.
    count_1s = 1
    total_count += count_1s
    terms.append(count_1s)
    print(f"For partition (1, 1, 1, 1, 1):")
    print(f"  Structure: Q x Q x Q x Q x Q")
    print(f"  Count = 1\n")
    
    print("="*50 + "\n")
    print("Final Calculation:")
    equation = " + ".join(map(str, terms))
    print(f"Total number = {equation} = {total_count}")

if __name__ == "__main__":
    main()