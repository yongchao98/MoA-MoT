import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Given parameters from the problem description
    b_p = 80  # Number of pages in relation P
    b_q = 65  # Number of pages in relation Q
    m = 15    # Number of available memory buffer pages

    # For BNLJ, M-2 pages are used for the outer relation's block
    buffer_for_outer = m - 2

    # --- Case 1: P is the outer relation ---
    # Number of times we need to loop over the inner relation Q
    num_blocks_p = math.ceil(b_p / buffer_for_outer)
    # Total I/O cost = cost to read P + (number of P blocks * cost to read Q)
    cost_p_outer = b_p + (num_blocks_p * b_q)

    # --- Case 2: Q is the outer relation ---
    # Number of times we need to loop over the inner relation P
    num_blocks_q = math.ceil(b_q / buffer_for_outer)
    # Total I/O cost = cost to read Q + (number of Q blocks * cost to read P)
    cost_q_outer = b_q + (num_blocks_q * b_p)
    
    print("### Block Nested Loop Join Cost Analysis ###")
    print(f"Pages in P, B(P) = {b_p}")
    print(f"Pages in Q, B(Q) = {b_q}")
    print(f"Memory buffer pages, M = {m}")
    print(f"Buffer pages for outer relation block = M - 2 = {buffer_for_outer}\n")
    
    print("--- Calculating Cost with P as Outer Relation ---")
    print(f"Number of blocks for P = ceil({b_p} / {buffer_for_outer}) = {num_blocks_p}")
    print(f"Cost = B(P) + (num_blocks_p * B(Q)) = {b_p} + {num_blocks_p} * {b_q} = {cost_p_outer}\n")

    print("--- Calculating Cost with Q as Outer Relation ---")
    print(f"Number of blocks for Q = ceil({b_q} / {buffer_for_outer}) = {num_blocks_q}")
    print(f"Cost = B(Q) + (num_blocks_q * B(P)) = {b_q} + {num_blocks_q} * {b_p} = {cost_q_outer}\n")

    # Determine the minimum cost and print the final result
    print("--- Conclusion ---")
    if cost_p_outer < cost_q_outer:
        print("Minimum cost is achieved with P as the outer relation.")
        print(f"Minimum Cost Equation: {b_p} + {num_blocks_p} * {b_q} = {cost_p_outer}")
        final_answer = cost_p_outer
    else:
        print("Minimum cost is achieved with Q as the outer relation.")
        print(f"Minimum Cost Equation: {b_q} + {num_blocks_q} * {b_p} = {cost_q_outer}")
        final_answer = cost_q_outer
    
    # Do not remove the following line
    print(f"\n<<<{final_answer}>>>")

calculate_bnlj_cost()