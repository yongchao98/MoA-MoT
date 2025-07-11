import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def main():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.
    """
    # T(k) is the number of non-isomorphic connected functional graphs on k vertices.
    # These values are pre-calculated.
    T = {1: 1, 2: 2, 3: 4, 4: 9}

    print("The total number of equivalence classes is the sum of counts for each partition of 4.")
    print("-" * 70)

    # Partition {4}
    # One component of size 4.
    count_4 = T[4]
    print(f"For partition {{4}}: there are T(4) = {count_4} classes.")

    # Partition {3, 1}
    # One component of size 3 and one of size 1.
    count_3_1 = T[3] * T[1]
    print(f"For partition {{3, 1}}: there are T(3) * T(1) = {T[3]} * {T[1]} = {count_3_1} classes.")

    # Partition {2, 2}
    # Two components of size 2. Choose 2 from T(2) types with replacement.
    count_2_2 = combinations(T[2] + 2 - 1, 2)
    print(f"For partition {{2, 2}}: there are C(T(2)+2-1, 2) = C({T[2]}+1, 2) = C(3, 2) = {count_2_2} classes.")

    # Partition {2, 1, 1}
    # One of size 2, two of size 1.
    count_2_1_1 = combinations(T[2] + 1 - 1, 1) * combinations(T[1] + 2 - 1, 2)
    print(f"For partition {{2, 1, 1}}: there are C(T(2)+1-1, 1) * C(T(1)+2-1, 2) = {T[2]} * 1 = {count_2_1_1} classes.")

    # Partition {1, 1, 1, 1}
    # Four components of size 1.
    count_1_1_1_1 = combinations(T[1] + 4 - 1, 4)
    print(f"For partition {{1, 1, 1, 1}}: there are C(T(1)+4-1, 4) = C({T[1]}+3, 4) = C(4, 4) = {count_1_1_1_1} class.")
    
    print("-" * 70)

    total_classes = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1
    
    print(f"Total = {count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_classes}")
    print(f"\nSo, there are {total_classes} elements of E represented by F.")

if __name__ == "__main__":
    main()