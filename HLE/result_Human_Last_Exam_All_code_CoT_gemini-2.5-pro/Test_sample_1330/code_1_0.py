import math

def solve():
    """
    Calculates the minimum cost of a Block Nested Loop Join (BNLJ)
    and prints the step-by-step calculation for the optimal case.
    """
    # Given parameters
    bp = 80  # Pages in relation P
    bq = 65  # Pages in relation Q
    m = 15   # Memory buffer pages

    # The buffer pages available for the outer relation block
    buffer_for_outer = m - 2

    # --- Scenario 1: P is the outer relation ---
    # Number of times we need to read the inner relation Q completely
    chunks_p = math.ceil(bp / buffer_for_outer)
    # Total cost = cost to read P once + (number of P chunks * cost to read Q)
    cost_p_outer = bp + chunks_p * bq

    # --- Scenario 2: Q is the outer relation ---
    # Number of times we need to read the inner relation P completely
    chunks_q = math.ceil(bq / buffer_for_outer)
    # Total cost = cost to read Q once + (number of Q chunks * cost to read P)
    cost_q_outer = bq + chunks_q * bp

    print("To find the minimum join cost, we evaluate two scenarios.\n")
    print(f"Given: B(P) = {bp}, B(Q) = {bq}, Memory pages M = {m}")
    print(f"Buffer pages for outer relation block = M - 2 = {m} - 2 = {buffer_for_outer}\n")

    print("Scenario 1: P as outer relation, Q as inner relation.")
    print(f"Cost = B(P) + ceil(B(P) / (M - 2)) * B(Q)")
    print(f"Cost = {bp} + ceil({bp} / {buffer_for_outer}) * {bq}")
    print(f"Cost = {bp} + {chunks_p} * {bq} = {cost_p_outer}\n")

    print("Scenario 2: Q as outer relation, P as inner relation.")
    print(f"Cost = B(Q) + ceil(B(Q) / (M - 2)) * B(P)")
    print(f"Cost = {bq} + ceil({bq} / {buffer_for_outer}) * {bp}")
    print(f"Cost = {bq} + {chunks_q} * {bp} = {cost_q_outer}\n")

    if cost_p_outer < cost_q_outer:
        min_cost = cost_p_outer
        outer_relation = "P"
        outer_pages = bp
        inner_pages = bq
        num_chunks = chunks_p
    else:
        min_cost = cost_q_outer
        outer_relation = "Q"
        outer_pages = bq
        inner_pages = bp
        num_chunks = chunks_q

    print(f"Comparing the two scenarios, the minimum cost is {min_cost} I/O operations.")
    print(f"This is achieved by using the smaller relation, {outer_relation}, as the outer relation.\n")
    print("The final calculation is:")
    print(f"{outer_pages} + {num_chunks} * {inner_pages} = {min_cost}")

    print(f"<<<{min_cost}>>>")

solve()