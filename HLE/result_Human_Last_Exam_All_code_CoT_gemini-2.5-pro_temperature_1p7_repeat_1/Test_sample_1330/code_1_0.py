import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Step 1: Define the given parameters.
    b_p = 80  # Pages in relation P
    b_q = 65  # Pages in relation Q
    m = 15    # Memory buffer pages

    # The BNLJ algorithm uses M-1 pages for the outer relation's block
    # and reserves one page for scanning the inner relation.
    buffer_for_outer = m - 1

    # Step 2: Calculate the cost with P as the outer relation.
    # The number of outer loops is the ceiling of B(P) / (M-1).
    num_p_blocks = math.ceil(b_p / buffer_for_outer)
    cost_p_outer = b_p + (num_p_blocks * b_q)

    # Step 3: Calculate the cost with Q as the outer relation.
    # The number of outer loops is the ceiling of B(Q) / (M-1).
    num_q_blocks = math.ceil(b_q / buffer_for_outer)
    cost_q_outer = b_q + (num_q_blocks * b_p)

    # Step 4: Determine the minimum cost and the optimal strategy.
    print("Calculating cost with P as the outer relation:")
    print(f"Cost = B(P) + (ceil(B(P) / (M - 1)) * B(Q))")
    print(f"Cost = {b_p} + (ceil({b_p} / ({m} - 1)) * {b_q})")
    print(f"Cost = {b_p} + ({num_p_blocks} * {b_q})")
    print(f"Cost = {b_p} + {num_p_blocks * b_q}")
    print(f"Cost = {cost_p_outer}\n")

    print("Calculating cost with Q as the outer relation:")
    print(f"Cost = B(Q) + (ceil(B(Q) / (M - 1)) * B(P))")
    print(f"Cost = {b_q} + (ceil({b_q} / ({m} - 1)) * {b_p})")
    print(f"Cost = {b_q} + ({num_q_blocks} * {b_p})")
    print(f"Cost = {b_q} + {num_q_blocks * b_p}")
    print(f"Cost = {cost_q_outer}\n")
    
    # The minimum cost is found by choosing the smaller relation (Q) as the outer relation.
    min_cost = min(cost_p_outer, cost_q_outer)
    
    print(f"The minimum I/O cost is achieved when the smaller relation, Q, is used as the outer relation.")
    print(f"Minimum Cost = {cost_q_outer}")

if __name__ == "__main__":
    calculate_bnlj_cost()