import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Given parameters
    b_p = 80  # Number of pages for relation P
    b_q = 65  # Number of pages for relation Q
    m = 15    # Number of available memory pages

    # The number of buffer pages available for the outer relation block
    # M-2: one page for the inner relation, one for the output buffer
    buffer_for_outer = m - 2

    # Case 1: P is the outer relation, Q is the inner relation
    num_chunks_p = math.ceil(b_p / buffer_for_outer)
    cost_p_outer = b_p + (num_chunks_p * b_q)

    # Case 2: Q is the outer relation, P is the inner relation
    num_chunks_q = math.ceil(b_q / buffer_for_outer)
    cost_q_outer = b_q + (num_chunks_q * b_p)
    
    print("--- Block Nested Loop Join Cost Calculation ---")
    print(f"Pages for P, B(P) = {b_p}")
    print(f"Pages for Q, B(Q) = {b_q}")
    print(f"Memory pages, M = {m}")
    print(f"Buffer pages for outer relation chunk = M - 2 = {buffer_for_outer}\n")

    print("Case 1: P as outer relation")
    print(f"Number of chunks of P = ceil(B(P) / (M-2)) = ceil({b_p} / {buffer_for_outer}) = {num_chunks_p}")
    print(f"Cost = B(P) + (Chunks of P * B(Q)) = {b_p} + ({num_chunks_p} * {b_q}) = {cost_p_outer}\n")
    
    print("Case 2: Q as outer relation")
    print(f"Number of chunks of Q = ceil(B(Q) / (M-2)) = ceil({b_q} / {buffer_for_outer}) = {num_chunks_q}")
    print(f"Cost = B(Q) + (Chunks of Q * B(P)) = {b_q} + ({num_chunks_q} * {b_p}) = {cost_q_outer}\n")
    
    print("--- Conclusion ---")
    if cost_q_outer < cost_p_outer:
        print("The minimum cost is achieved when Q is the outer relation.")
        print("Final minimum cost equation:")
        print(f"{b_q} + {num_chunks_q} * {b_p} = {cost_q_outer}")
    else:
        print("The minimum cost is achieved when P is the outer relation.")
        print("Final minimum cost equation:")
        print(f"{b_p} + {num_chunks_p} * {b_q} = {cost_p_outer}")

calculate_bnlj_cost()
<<<465>>>