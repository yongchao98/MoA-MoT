def find_local_flow_solution():
    """
    This function checks for the smallest integer k > 1 for which the local
    flow conservation law, a + b + c = 0, can be satisfied at a vertex
    in a 3-regular graph. The flow values a, b, c must be in the set
    {+/-1, ..., +/-(k-1)}.
    """

    # We start checking from k=2, as k=1 results in an empty set of values.
    for k in range(2, 6):
        print(f"--- Checking for k={k} ---")
        flow_values = list(range(-(k - 1), 0)) + list(range(1, k))
        
        solution_found = False
        # Iterate through all combinations of three flow values
        for a in flow_values:
            for b in flow_values:
                for c in flow_values:
                    if a + b + c == 0:
                        print(f"A local solution exists for k={k}.")
                        print("For example, the following equation using these values sums to zero:")
                        # As requested, output each number in the final equation.
                        print(f"{a} + {b} + {c} = 0")
                        solution_found = True
                        break  # Exit inner loop
                if solution_found:
                    break  # Exit middle loop
            if solution_found:
                break  # Exit outer loop
        
        if not solution_found:
            print(f"No local solution found for k={k}.")
        
        # We only need to show the first k that works locally.
        if solution_found:
            print("\nWhile a local solution exists for k=3, this does not guarantee a global solution for the entire graph.")
            print("Graph theory theorems show that for certain graphs, a k of 5 is required.")
            break

find_local_flow_solution()