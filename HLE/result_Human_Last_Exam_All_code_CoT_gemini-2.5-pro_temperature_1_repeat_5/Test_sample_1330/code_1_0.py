import math

def solve():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join (BNLJ).
    """
    # Given parameters
    b_p = 80  # Number of pages for relation P
    b_q = 65  # Number of pages for relation Q
    m = 15    # Number of available memory buffer pages

    # Buffer pages available for the outer relation
    # M-2 because one page is for the inner relation and one is for the output buffer.
    buffer_for_outer = m - 2

    print("--- Calculating I/O Cost for Block Nested Loop Join ---")
    print(f"Pages for P, B(P) = {b_p}")
    print(f"Pages for Q, B(Q) = {b_q}")
    print(f"Memory Buffer Pages, M = {m}")
    print(f"Buffer pages available for the outer relation blocks = M - 2 = {buffer_for_outer}\n")

    # Case 1: P is the outer relation, Q is the inner relation
    # Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))
    num_p_chunks = math.ceil(b_p / buffer_for_outer)
    cost_p_outer = b_p + num_p_chunks * b_q

    print("--- Case 1: P is the outer relation, Q is the inner ---")
    print(f"Number of chunks for P = ceil(B(P) / (M - 2)) = ceil({b_p} / {buffer_for_outer}) = {num_p_chunks}")
    print(f"Cost = B(P) + (Number of P chunks * B(Q))")
    print(f"Cost = {b_p} + ({num_p_chunks} * {b_q}) = {b_p} + {num_p_chunks * b_q} = {cost_p_outer}\n")

    # Case 2: Q is the outer relation, P is the inner relation
    # Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))
    num_q_chunks = math.ceil(b_q / buffer_for_outer)
    cost_q_outer = b_q + num_q_chunks * b_p

    print("--- Case 2: Q is the outer relation, P is the inner ---")
    print(f"Number of chunks for Q = ceil(B(Q) / (M - 2)) = ceil({b_q} / {buffer_for_outer}) = {num_q_chunks}")
    print(f"Cost = B(Q) + (Number of Q chunks * B(P))")
    print(f"Cost = {b_q} + ({num_q_chunks} * {b_p}) = {b_q} + {num_q_chunks * b_p} = {cost_q_outer}\n")

    # Determine the minimum cost
    min_cost = min(cost_p_outer, cost_q_outer)

    print("--- Conclusion ---")
    print(f"Comparing the two costs: Cost(P outer) = {cost_p_outer}, Cost(Q outer) = {cost_q_outer}")
    print(f"The minimum cost is achieved by making the smaller relation (Q) the outer relation.")
    print(f"Minimum I/O Cost = {min_cost}")
    print(f"<<<{min_cost}>>>")

solve()