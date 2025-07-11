import math

def calculate_bnlj_cost():
    """
    Calculates the minimum cost for a Block Nested Loop Join (BNLJ) given two relations and a buffer size.
    """
    # Given parameters
    b_p = 80  # Number of pages for relation P
    b_q = 65  # Number of pages for relation Q
    m = 15    # Number of available memory buffer pages

    # The number of buffer pages available to hold blocks of the outer relation.
    # We reserve one page for the inner relation's input and one for the output buffer.
    buffer_for_outer = m - 2

    print("This script calculates the minimum I/O cost for a Block Nested Loop Join.")
    print("The filtering conditions do not reduce I/O costs as the join is performed without materialization.")
    print("-" * 40)
    print(f"Given parameters:")
    print(f"  Pages in P, B(P) = {b_p}")
    print(f"  Pages in Q, B(Q) = {b_q}")
    print(f"  Memory pages, M = {m}")
    print(f"  Buffer pages for outer relation = M - 2 = {m} - 2 = {buffer_for_outer}")
    print("-" * 40)

    # --- Case 1: P as outer relation, Q as inner relation ---
    # Number of chunks P is divided into
    num_chunks_p = math.ceil(b_p / buffer_for_outer)
    # Cost = B(P) + (Number of P chunks) * B(Q)
    cost_p_outer = b_p + num_chunks_p * b_q

    print("Scenario 1: P is the outer relation.")
    print(f"Number of outer relation chunks = ceil(B(P) / (M - 2)) = ceil({b_p} / {buffer_for_outer}) = {num_chunks_p}")
    print("Cost formula: B(P) + (Number of chunks of P) * B(Q)")
    print(f"Cost = {b_p} + {num_chunks_p} * {b_q} = {cost_p_outer}")
    print("-" * 40)

    # --- Case 2: Q as outer relation, P as inner relation ---
    # Number of chunks Q is divided into
    num_chunks_q = math.ceil(b_q / buffer_for_outer)
    # Cost = B(Q) + (Number of Q chunks) * B(P)
    cost_q_outer = b_q + num_chunks_q * b_p
    
    print("Scenario 2: Q is the outer relation.")
    print(f"Number of outer relation chunks = ceil(B(Q) / (M - 2)) = ceil({b_q} / {buffer_for_outer}) = {num_chunks_q}")
    print("Cost formula: B(Q) + (Number of chunks of Q) * B(P)")
    print(f"Cost = {b_q} + {num_chunks_q} * {b_p} = {cost_q_outer}")
    print("-" * 40)

    # Determine the minimum cost
    min_cost = min(cost_p_outer, cost_q_outer)
    
    print(f"To minimize I/O operations, we choose the smaller relation (Q) as the outer relation.")
    print(f"Minimum Cost = min({cost_p_outer}, {cost_q_outer}) = {min_cost} I/O operations.")

calculate_bnlj_cost()