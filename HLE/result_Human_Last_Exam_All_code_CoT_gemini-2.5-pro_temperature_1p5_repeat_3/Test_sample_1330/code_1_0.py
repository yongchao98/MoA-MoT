import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Given parameters
    b_p = 80  # Number of pages for relation P
    b_q = 65  # Number of pages for relation Q
    m = 15    # Number of available memory buffer pages

    # The block nested loop join algorithm uses M-1 pages for the outer relation block
    # and 1 page for the inner relation scan.

    # Scenario 1: P is the outer relation, Q is the inner relation
    num_blocks_p = math.ceil(b_p / (m - 1))
    cost_p_outer = b_p + (num_blocks_p * b_q)

    # Scenario 2: Q is the outer relation, P is the inner relation
    num_blocks_q = math.ceil(b_q / (m - 1))
    cost_q_outer = b_q + (num_blocks_q * b_p)

    # Determine the minimum cost and the corresponding calculation
    if cost_p_outer < cost_q_outer:
        min_cost = cost_p_outer
        outer_relation_pages = b_p
        inner_relation_pages = b_q
        num_outer_blocks = num_blocks_p
        outer_relation_name = "P"
        inner_relation_name = "Q"
    else:
        min_cost = cost_q_outer
        outer_relation_pages = b_q
        inner_relation_pages = b_p
        num_outer_blocks = num_blocks_q
        outer_relation_name = "Q"
        inner_relation_name = "P"
        
    print(f"To minimize I/O cost, we choose the smaller relation as the outer relation.")
    print(f"The minimum cost is achieved when {outer_relation_name} is the outer relation.")
    print("\n--- Calculation for Minimum Cost ---")
    print(f"Cost = B({outer_relation_name}) + (ceil(B({outer_relation_name}) / (M - 1)) * B({inner_relation_name}))")
    print(f"Cost = {outer_relation_pages} + (ceil({outer_relation_pages} / ({m} - 1)) * {inner_relation_pages})")
    print(f"Cost = {outer_relation_pages} + ({num_outer_blocks} * {inner_relation_pages})")
    print(f"Cost = {outer_relation_pages} + {num_outer_blocks * inner_relation_pages}")
    print(f"Cost = {min_cost}")

calculate_bnlj_cost()
<<<465>>>