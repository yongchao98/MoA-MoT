import math

def solve_bnlj_cost():
    """
    Calculates the minimum cost for a Block Nested Loop Join given relation sizes
    and memory buffer size. It then prints the step-by-step calculation for the
    optimal join order.
    """
    # Given parameters
    b_p = 80  # Number of pages for relation P
    b_q = 65  # Number of pages for relation Q
    m = 15    # Number of available memory buffer pages

    # In BNLJ, we use M-2 pages for the outer relation block, 1 for the inner,
    # and 1 for the output buffer.
    buffer_for_outer = m - 2

    # --- Scenario 1: P as outer relation, Q as inner ---
    num_outer_loops_p = math.ceil(b_p / buffer_for_outer)
    cost_p_outer = b_p + (num_outer_loops_p * b_q)

    # --- Scenario 2: Q as outer relation, P as inner ---
    num_outer_loops_q = math.ceil(b_q / buffer_for_outer)
    cost_q_outer = b_q + (num_outer_loops_q * b_p)

    # Determine the minimum cost and print the detailed calculation for the optimal scenario.
    print("To minimize I/O cost, the smaller relation should be the outer relation.")
    
    if cost_p_outer <= cost_q_outer:
        final_cost = cost_p_outer
        print("Optimal Plan: Use P as the outer relation and Q as the inner relation.\n")
        print("Calculation using BNLJ cost formula: B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
        print(f"Cost = {b_p} + (ceil({b_p} / ({m} - 2)) * {b_q})")
        print(f"Cost = {b_p} + (ceil({b_p / buffer_for_outer:.2f}) * {b_q})")
        print(f"Cost = {b_p} + ({num_outer_loops_p} * {b_q})")
        print(f"Cost = {b_p} + {num_outer_loops_p * b_q}")
        print(f"Total minimum I/O cost = {final_cost}")
    else:
        final_cost = cost_q_outer
        print("Optimal Plan: Use Q (65 pages) as the outer relation and P (80 pages) as the inner relation.\n")
        print("Calculation using BNLJ cost formula: B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
        print(f"Cost = {b_q} + (ceil({b_q} / ({m} - 2)) * {b_p})")
        print(f"Cost = {b_q} + (ceil({b_q / buffer_for_outer:.2f}) * {b_p})")
        print(f"Cost = {b_q} + ({num_outer_loops_q} * {b_p})")
        print(f"Cost = {b_q} + {num_outer_loops_q * b_p}")
        print(f"Total minimum I/O cost = {final_cost}")

solve_bnlj_cost()
<<<465>>>