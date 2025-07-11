import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Given parameters
    b_p = 80  # Pages in relation P
    b_q = 65  # Pages in relation Q
    m = 15    # Memory buffer pages

    # In BNLJ, we use M-2 pages for the outer relation's block,
    # 1 page for the inner relation, and 1 page for the output buffer.
    buffer_for_outer = m - 2

    # --- Case 1: P is the outer relation ---
    num_outer_chunks_p = math.ceil(b_p / buffer_for_outer)
    cost_p_outer = b_p + num_outer_chunks_p * b_q

    # --- Case 2: Q is the outer relation ---
    num_outer_chunks_q = math.ceil(b_q / buffer_for_outer)
    cost_q_outer = b_q + num_outer_chunks_q * b_p

    print("--- Block Nested Loop Join Cost Calculation ---")
    print(f"Relation P pages B(P) = {b_p}")
    print(f"Relation Q pages B(Q) = {b_q}")
    print(f"Memory buffer pages M = {m}")
    print(f"Buffer pages for outer relation block = M - 2 = {buffer_for_outer}\n")

    print("Case 1: P as the outer relation")
    print(f"Cost = B(P) + (ceil(B(P) / (M - 2))) * B(Q)")
    print(f"Cost = {b_p} + (ceil({b_p} / {buffer_for_outer})) * {b_q}")
    print(f"Cost = {b_p} + {num_outer_chunks_p} * {b_q}")
    print(f"Cost = {b_p} + {num_outer_chunks_p * b_q}")
    print(f"Cost = {cost_p_outer} I/Os\n")

    print("Case 2: Q as the outer relation")
    print(f"Cost = B(Q) + (ceil(B(Q) / (M - 2))) * B(P)")
    print(f"Cost = {b_q} + (ceil({b_q} / {buffer_for_outer})) * {b_p}")
    print(f"Cost = {b_q} + {num_outer_chunks_q} * {b_p}")
    print(f"Cost = {b_q} + {num_outer_chunks_q * b_p}")
    print(f"Cost = {cost_q_outer} I/Os\n")

    print("--- Conclusion ---")
    print("To minimize I/O cost, we choose the smaller relation (Q) as the outer relation.")
    
    min_cost = min(cost_p_outer, cost_q_outer)
    
    print(f"Minimum Cost = {min_cost} I/Os")

if __name__ == '__main__':
    calculate_bnlj_cost()
    # The final answer is the minimum cost calculated
    print("\n<<<465>>>")
