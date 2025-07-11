import math

def combinations(n, k):
    """Calculates the binomial coefficient 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def count_structures():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.
    This corresponds to the number of non-isomorphic functional graphs on 4 vertices.
    """
    # c_k is the number of connected functional graphs on k vertices.
    # OEIS A000201
    c = {
        1: 1,
        2: 2,
        3: 4,
        4: 9
    }

    total_classes = 0

    # Case 1: Partition {4} -> one component of size 4
    # Choose 1 component from c_4 types.
    count1 = c[4]
    total_classes += count1
    print(f"Partition {{4}}: {count1} structures")

    # Case 2: Partition {3, 1} -> one of size 3, one of size 1
    # Choose 1 from c_3 and 1 from c_1.
    count2 = c[3] * c[1]
    total_classes += count2
    print(f"Partition {{3, 1}}: c(3) * c(1) = {c[3]} * {c[1]} = {count2} structures")

    # Case 3: Partition {2, 2} -> two components of size 2
    # Choose a multiset of 2 from c_2 types.
    count3 = combinations(c[2] + 2 - 1, 2)
    total_classes += count3
    print(f"Partition {{2, 2}}: C({c[2]}+2-1, 2) = {count3} structures")

    # Case 4: Partition {2, 1, 1} -> one of size 2, two of size 1
    # Choose 1 from c_2 types and a multiset of 2 from c_1 type.
    count4 = combinations(c[2], 1) * combinations(c[1] + 2 - 1, 2)
    total_classes += count4
    print(f"Partition {{2, 1, 1}}: C({c[2]}, 1) * C({c[1]}+2-1, 2) = {combinations(c[2], 1)} * {combinations(c[1] + 2 - 1, 2)} = {count4} structures")

    # Case 5: Partition {1, 1, 1, 1} -> four components of size 1
    # Choose a multiset of 4 from c_1 type.
    count5 = combinations(c[1] + 4 - 1, 4)
    total_classes += count5
    print(f"Partition {{1, 1, 1, 1}}: C({c[1]}+4-1, 4) = {count5} structure")
    
    print("-" * 20)
    print(f"Total number of classes = {count1} + {count2} + {count3} + {count4} + {count5}")
    print(f"Final Answer: {total_classes}")

count_structures()