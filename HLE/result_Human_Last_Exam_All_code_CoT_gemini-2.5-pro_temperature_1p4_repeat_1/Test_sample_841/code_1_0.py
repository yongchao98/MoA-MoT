import math

def solve():
    """
    Solves the problem by counting the number of HFDs in the specified set of rings.
    """

    # Part 1: Count maximal orders that are HFDs.
    # An order O_K is an HFD iff its class number h(-d) <= 2.
    # We use the known lists of such d's.

    # d values for which h(-d) = 1 (Heegner numbers)
    h1_d = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    num_h1 = len(h1_d)

    # d values for which h(-d) = 2
    h2_d = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]
    num_h2 = len(h2_d)

    num_maximal_hfd = num_h1 + num_h2
    print(f"Number of maximal orders (rings of integers of Q(sqrt(-d))) that are HFDs:")
    print(f"Found {num_h1} rings with class number 1 and {num_h2} rings with class number 2.")
    print(f"Total is {num_h1} + {num_h2} = {num_maximal_hfd}.")
    print("-" * 20)

    # Part 2: Count non-maximal orders Z[sqrt(-d)] that are HFDs.
    # This applies only when d = 3 (mod 4). The order is an HFD if its class number, h(R), is <= 2.
    # We only need to check d values where the maximal order class number h(-d) is <= 2.
    
    non_maximal_hfd_d = []
    
    # Combine the lists of d's where h(-d) <= 2.
    candidates_d = h1_d + h2_d

    for d in candidates_d:
        # We are only interested in the rings Z[sqrt(-d)] that are not integrally closed,
        # which means d must be 3 (mod 4).
        if d % 4 == 3:
            # Determine h_K (class number of the maximal order)
            h_K = 1 if d in h1_d else 2
            
            # The class number h(R) of the non-maximal order R=Z[sqrt(-d)] is calculated.
            # A known result states:
            # For d=3, h(R) = 1.
            # For d > 3, d = 3 (mod 8), h(R) = 3*h_K.
            # For d > 3, d = 7 (mod 8), h(R) = h_K.
            
            h_R = 0
            if d == 3:
                h_R = 1
            elif d > 3 and d % 8 == 3:
                h_R = 3 * h_K
            elif d % 8 == 7: # d > 3 is implied, as 3 % 8 is not 7.
                h_R = h_K
            
            if h_R <= 2:
                non_maximal_hfd_d.append(d)

    num_non_maximal_hfd = len(non_maximal_hfd_d)

    print(f"Number of non-maximal orders (Z[sqrt(-d)] with d=3(mod 4)) that are HFDs:")
    print(f"Found {num_non_maximal_hfd} such rings, for d in {non_maximal_hfd_d}.")
    print("-" * 20)
    
    # Part 3: Calculate the total size of the subset.
    # The two sets of rings (maximal orders and the specified non-maximal orders) are disjoint.
    total_hfd_rings = num_maximal_hfd + num_non_maximal_hfd

    print("The total size of the subset of HFDs is the sum of these two counts.")
    print("Final Equation:")
    print(f"{num_maximal_hfd} + {num_non_maximal_hfd} = {total_hfd_rings}")
    
    return total_hfd_rings

# Run the solver and capture the final answer
final_answer = solve()
print(f"\n<<< {final_answer} >>>")
