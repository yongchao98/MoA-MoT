def solve_and_print():
    """
    This function calculates the size of the specified set of number rings
    for which prime factorizations have unique lengths.
    """
    
    # Step 1: Define the square-free integers d > 0 for which the class number h
    # of the maximal order O(Q(sqrt(-d))) is 1 or 2. These are known results.
    d_h1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    d_h2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    # Combine the lists for all d where the maximal order is an HFD.
    all_d_maximal_hfd = sorted(d_h1 + d_h2)

    # Step 2: Count the maximal orders with class number <= 2.
    # These correspond to every d in the lists above.
    num_maximal_hfd = len(all_d_maximal_hfd)

    print("Part 1: Counting Maximal Orders with Unique Factorization Lengths (HFDs)")
    print("These are the rings of integers O(Q(sqrt(-d))) with class number h(O_d) <= 2.")
    print(f"There are {len(d_h1)} values of d for which h(O_d) = 1.")
    print(f"There are {len(d_h2)} values of d for which h(O_d) = 2.")
    print(f"Total number of such maximal orders is {len(d_h1)} + {len(d_h2)} = {num_maximal_hfd}.\n")

    # Step 3: Count non-maximal orders Z[sqrt(-d)] with class number <= 2.
    # These rings are non-maximal only when d = 3 (mod 4).
    # We find which of these have a class number <= 2.
    
    print("Part 2: Counting Non-Maximal Orders with Unique Factorization Lengths (HFDs)")
    print("These are rings Z[sqrt(-d)] where d=3(mod 4) and class number h(R_d) <= 2.")
    
    non_maximal_hfd_d = []
    
    # We only need to check d from our known lists, as h(R_d) is related to h(O_d).
    d_candidates = [d for d in all_d_maximal_hfd if d % 4 == 3]

    for d in d_candidates:
        h_Od = 1 if d in d_h1 else 2
        h_Rd = 0
        
        # We use the known formula for the class number h(R_d) of a non-maximal order.
        if d == 3:
            # Special case for d=3 due to having more units.
            # h(R_3) = (h(O_3)/3) * 2 * (1 - (-3|2)/2), where (-3|2) = -1
            h_Rd = 1
        else: # d > 3 and d = 3 (mod 4)
            if d % 8 == 3:
                # h(R_d) = h(O_d) * 3
                h_Rd = 3 * h_Od
            elif d % 8 == 7:
                # h(R_d) = h(O_d)
                h_Rd = h_Od

        if h_Rd <= 2:
            non_maximal_hfd_d.append(d)

    num_non_maximal_hfd = len(non_maximal_hfd_d)
    print(f"The values of d for which Z[sqrt(-d)] is a non-maximal HFD are: {non_maximal_hfd_d}.")
    print(f"Total number of such non-maximal orders is {num_non_maximal_hfd}.\n")

    # Step 4: Calculate the total size of the union.
    # The sets of maximal and non-maximal rings are disjoint, so we add their sizes.
    total_size = num_maximal_hfd + num_non_maximal_hfd

    print("Final Calculation")
    print("The size of the union is the sum of the counts for the two disjoint sets.")
    print(f"Total Size = (Number of maximal HFDs) + (Number of non-maximal HFDs)")
    print(f"Total Size = {num_maximal_hfd} + {num_non_maximal_hfd} = {total_size}")


if __name__ == '__main__':
    solve_and_print()