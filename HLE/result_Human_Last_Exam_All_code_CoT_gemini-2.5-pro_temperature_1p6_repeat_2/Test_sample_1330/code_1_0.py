import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Given parameters
    b_p = 80  # Pages in relation P
    b_q = 65  # Pages in relation Q
    m = 15    # Available memory buffer pages

    # Reserved buffer pages (1 for inner scan, 1 for output)
    reserved_pages = 2
    available_buffer_for_outer = m - reserved_pages

    # --- Scenario 1: P is the outer relation ---
    num_chunks_p_outer = math.ceil(b_p / available_buffer_for_outer)
    cost_p_outer = b_p + num_chunks_p_outer * b_q

    # --- Scenario 2: Q is the outer relation ---
    num_chunks_q_outer = math.ceil(b_q / available_buffer_for_outer)
    cost_q_outer = b_q + num_chunks_q_outer * b_p

    # Determine the minimum cost and print the calculation steps for the optimal scenario
    if cost_q_outer < cost_p_outer:
        print("The minimum cost is achieved by using Q as the outer relation and P as the inner relation.")
        print("The BNLJ cost formula is: B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
        print(f"Calculation:")
        print(f"Cost = {b_q} + (ceil({b_q} / ({m} - {reserved_pages})) * {b_p})")
        print(f"Cost = {b_q} + (ceil({b_q} / {available_buffer_for_outer}) * {b_p})")
        # intermediate step for ceil calculation
        ceil_val = b_q / available_buffer_for_outer
        print(f"Cost = {b_q} + (ceil({ceil_val}) * {b_p})")
        print(f"Cost = {b_q} + ({num_chunks_q_outer} * {b_p})")
        inner_scan_cost = num_chunks_q_outer * b_p
        print(f"Cost = {b_q} + {inner_scan_cost}")
        final_cost = b_q + inner_scan_cost
        print(f"Total Minimum Cost = {final_cost} I/O operations")
        
    else:
        print("The minimum cost is achieved by using P as the outer relation and Q as the inner relation.")
        print("The BNLJ cost formula is: B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
        print(f"Calculation:")
        print(f"Cost = {b_p} + (ceil({b_p} / ({m} - {reserved_pages})) * {b_q})")
        print(f"Cost = {b_p} + (ceil({b_p} / {available_buffer_for_outer}) * {b_q})")
        # intermediate step for ceil calculation
        ceil_val = b_p / available_buffer_for_outer
        print(f"Cost = {b_p} + (ceil({ceil_val:.2f}) * {b_q})")
        print(f"Cost = {b_p} + ({num_chunks_p_outer} * {b_q})")
        inner_scan_cost = num_chunks_p_outer * b_q
        print(f"Cost = {b_p} + {inner_scan_cost}")
        final_cost = b_p + inner_scan_cost
        print(f"Total Minimum Cost = {final_cost} I/O operations")

calculate_bnlj_cost()
<<<465>>>