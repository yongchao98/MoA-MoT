def solve_unique_length_factorizations():
    """
    Calculates the size of the set of imaginary quadratic rings where prime factorizations
    have unique lengths.
    """

    # These are the known square-free d > 0 for which the class number h of Q(sqrt(-d))
    # is 1 or 2. For the ring of integers of Q(sqrt(-d)), being a Half-Factorial
    # Domain (HFD) is equivalent to h=1 or h=2.
    d_h1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    d_h2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    # The union of rings can be partitioned into three disjoint sets.
    # We count the number of HFDs in each set.

    # ---- Subset 1: Rings of integers Z[sqrt(-d)] for d = 1, 2 (mod 4) ----
    # These are HFDs if the class number is 1 or 2.
    h1_subset1 = [d for d in d_h1 if d % 4 == 1 or d % 4 == 2]
    h2_subset1 = [d for d in d_h2 if d % 4 == 1 or d % 4 == 2]
    count1 = len(h1_subset1) + len(h2_subset1)

    # ---- Subset 2: Rings of integers Z[(1+sqrt(-d))/2] for d = 3 (mod 4) ----
    # These are HFDs if the class number is 1 or 2.
    h1_subset2 = [d for d in d_h1 if d % 4 == 3]
    h2_subset2 = [d for d in d_h2 if d % 4 == 3]
    count2 = len(h1_subset2) + len(h2_subset2)

    # ---- Subset 3: Non-maximal orders Z[sqrt(-d)] for d = 3 (mod 4) ----
    # For these orders (with d square-free), they are HFDs only for d=3 and d=7.
    # We check they satisfy the condition d=3 (mod 4).
    d_subset3 = []
    if 3 % 4 == 3:
        d_subset3.append(3)
    if 7 % 4 == 3:
        d_subset3.append(7)
    count3 = len(d_subset3)

    # The total number is the sum of counts from the disjoint subsets.
    total_size = count1 + count2 + count3

    print("The problem asks for the size of a set of number rings that are Half-Factorial Domains (HFDs).")
    print("The set of rings can be partitioned into three disjoint subsets:\n")
    print(f"1. Rings of integers for d = 1, 2 (mod 4): The number of HFDs is {len(h1_subset1)} (from h=1) + {len(h2_subset1)} (from h=2) = {count1}.")
    print(f"2. Rings of integers for d = 3 (mod 4): The number of HFDs is {len(h1_subset2)} (from h=1) + {len(h2_subset2)} (from h=2) = {count2}.")
    print(f"3. Non-maximal orders for d = 3 (mod 4): The number of HFDs is {count3}.\n")

    print("The total size of the set of HFDs is the sum of these counts.")
    print(f"Total size = {count1} + {count2} + {count3} = {total_size}")


solve_unique_length_factorizations()