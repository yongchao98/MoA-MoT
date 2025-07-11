import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Step 1: Define the given parameters
    b_p = 80  # Pages in relation P
    b_q = 65  # Pages in relation Q
    m = 15    # Available memory pages

    print("This script calculates the minimum I/O cost for a Block Nested Loop Join.")
    print("The filter conditions are applied on-the-fly and do not reduce the I/O cost without materialization.")
    print("\n--- Given Parameters ---")
    print(f"Pages in P, B(P) = {b_p}")
    print(f"Pages in Q, B(Q) = {b_q}")
    print(f"Memory Buffer Pages, M = {m}")

    # For BNLJ, M-2 pages are used for the outer relation block
    buffer_for_outer = m - 2
    print(f"Buffer pages for the outer relation block = M - 2 = {m} - 2 = {buffer_for_outer}")

    # --- Scenario 1: P as the outer relation ---
    print("\n--- Scenario 1: P as Outer Relation ---")
    print("Cost = B(P) + (ceil(B(P) / (M-2)) * B(Q))")
    
    # Calculate the number of blocks for the outer relation P
    num_blocks_p = math.ceil(b_p / buffer_for_outer)
    cost_p_outer = b_p + num_blocks_p * b_q
    
    print(f"Number of blocks for P = ceil({b_p} / {buffer_for_outer}) = {num_blocks_p}")
    print(f"Cost = {b_p} + ({num_blocks_p} * {b_q}) = {b_p} + {num_blocks_p * b_q} = {cost_p_outer}")

    # --- Scenario 2: Q as the outer relation ---
    print("\n--- Scenario 2: Q as Outer Relation ---")
    print("Cost = B(Q) + (ceil(B(Q) / (M-2)) * B(P))")

    # Calculate the number of blocks for the outer relation Q
    num_blocks_q = math.ceil(b_q / buffer_for_outer)
    cost_q_outer = b_q + num_blocks_q * b_p
    
    print(f"Number of blocks for Q = ceil({b_q} / {buffer_for_outer}) = {num_blocks_q}")
    print(f"Cost = {b_q} + ({num_blocks_q} * {b_p}) = {b_q} + {num_blocks_q * b_p} = {cost_q_outer}")

    # --- Conclusion: Determine the minimum cost ---
    print("\n--- Conclusion ---")
    print(f"Comparing the two costs: {cost_p_outer} (P as outer) and {cost_q_outer} (Q as outer).")
    
    if cost_p_outer < cost_q_outer:
        min_cost = cost_p_outer
        print(f"The minimum cost is {min_cost} I/Os, achieved with P as the outer relation.")
        print("\nFinal minimum cost calculation:")
        print(f"{b_p} + ({num_blocks_p} * {b_q}) = {min_cost}")
    else:
        min_cost = cost_q_outer
        print(f"The minimum cost is {min_cost} I/Os, achieved with Q as the outer relation.")
        print("\nFinal minimum cost calculation:")
        print(f"{b_q} + ({num_blocks_q} * {b_p}) = {min_cost}")

    # The final answer in the required format
    # This part will not be printed to the user console, it's the final answer marker.
    return min_cost

# Execute the function to print the calculation steps
min_cost_value = calculate_bnlj_cost()
# print(f"\n<<<{min_cost_value}>>>") # This is for the final answer extraction