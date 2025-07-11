import math

def combinations(n, k):
    """Calculates the binomial coefficient 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def main():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction outside the prime 2.
    """
    # a_d is the number of number fields of degree d unramified outside the prime 2.
    # This data is based on results from algebraic number theory and databases like LMFDB.
    a = {
        1: 1,  # Q
        2: 3,  # Q(sqrt(-1)), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 9,
        5: 0
    }

    print("This program calculates the number of quintic étale algebras over Q unramified outside prime 2.")
    print("This count is equivalent to the number of isomorphism classes of del Pezzo surfaces of degree 5")
    print("over Q with good reduction everywhere except possibly at the prime 2.\n")
    print("The number of available number fields of degree d unramified outside 2 are:")
    for d, count in a.items():
        print(f"a_{d}: {count}")
    print("-" * 40)
    print("Calculating the number of algebras for each partition of 5:\n")
    
    total_count = 0
    
    # Partition 5: A single quintic field
    count_5 = a[5]
    total_count += count_5
    print(f"Partition [5]: A single field of degree 5.")
    print(f"Number of choices = a_5 = {a[5]}")
    print(f"Contribution = {count_5}\n")

    # Partition 4+1: A quartic field and a degree 1 field
    count_4_1 = a[4] * a[1]
    total_count += count_4_1
    print(f"Partition [4, 1]: A field of degree 4 and a field of degree 1.")
    print(f"Number of choices = a_4 * a_1 = {a[4]} * {a[1]} = {count_4_1}")
    print(f"Contribution = {count_4_1}\n")

    # Partition 3+2: A cubic field and a quadratic field
    count_3_2 = a[3] * a[2]
    total_count += count_3_2
    print(f"Partition [3, 2]: A field of degree 3 and a field of degree 2.")
    print(f"Number of choices = a_3 * a_2 = {a[3]} * {a[2]} = {count_3_2}")
    print(f"Contribution = {count_3_2}\n")

    # Partition 3+1+1: A cubic field and two degree 1 fields
    count_3_1_1 = a[3] * combinations(a[1] + 2 - 1, 2)
    total_count += count_3_1_1
    print(f"Partition [3, 1, 1]: A field of degree 3 and two fields of degree 1.")
    print(f"Number of choices = a_3 * C(a_1+2-1, 2) = {a[3]} * C({a[1]}+1, 2) = {count_3_1_1}")
    print(f"Contribution = {count_3_1_1}\n")

    # Partition 2+2+1: Two quadratic fields and a degree 1 field
    count_2_2_1 = combinations(a[2] + 2 - 1, 2) * a[1]
    total_count += count_2_2_1
    print(f"Partition [2, 2, 1]: Two fields of degree 2 and one field of degree 1.")
    print(f"The number of ways to choose 2 fields from {a[2]} with replacement is C({a[2]}+2-1, 2).")
    print(f"Number of choices = C(a_2+2-1, 2) * a_1 = C({a[2]}+1, 2) * {a[1]} = {combinations(a[2]+1, 2)} * {a[1]} = {count_2_2_1}")
    print(f"Contribution = {count_2_2_1}\n")

    # Partition 2+1+1+1: A quadratic field and three degree 1 fields
    count_2_1_1_1 = a[2] * combinations(a[1] + 3 - 1, 3)
    total_count += count_2_1_1_1
    print(f"Partition [2, 1, 1, 1]: One field of degree 2 and three fields of degree 1.")
    print(f"Number of choices = a_2 * C(a_1+3-1, 3) = {a[2]} * C({a[1]}+2, 3) = {a[2]} * {combinations(a[1]+2, 3)} = {count_2_1_1_1}")
    print(f"Contribution = {count_2_1_1_1}\n")

    # Partition 1+1+1+1+1: Five degree 1 fields
    count_1_1_1_1_1 = combinations(a[1] + 5 - 1, 5)
    total_count += count_1_1_1_1_1
    print(f"Partition [1, 1, 1, 1, 1]: Five fields of degree 1.")
    print(f"Number of choices = C(a_1+5-1, 5) = C({a[1]}+4, 5) = {combinations(a[1]+4, 5)} = {count_1_1_1_1_1}")
    print(f"Contribution = {count_1_1_1_1_1}\n")

    print("-" * 40)
    print(f"Total number of isomorphism classes = {count_5} + {count_4_1} + {count_3_2} + {count_3_1_1} + {count_2_2_1} + {count_2_1_1_1} + {count_1_1_1_1_1}")
    print(f"Total = {total_count}")
    
    # Final answer in the required format
    # The final answer is an integer, so no need for a complex format
    print("\n<<<19>>>")

if __name__ == "__main__":
    main()