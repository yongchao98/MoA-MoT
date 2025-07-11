import math

def calculate_bnlj_cost():
    """
    Calculates the minimum cost for a Block Nested Loop Join.

    The problem states:
    - Relation P pages B(P): 80
    - Relation Q pages B(Q): 65
    - Memory buffer pages M: 15
    - Selections must be applied on-the-fly without materialization.

    The cost formula for BNLJ is:
    Cost = B(outer) + (ceil(B(outer) / (M - 1)) * B(inner))
    
    We calculate the cost for two scenarios:
    1. P is the outer relation, Q is the inner relation.
    2. Q is the outer relation, P is the inner relation.
    
    The minimum of these two costs is the answer.
    """
    
    # Given parameters
    b_p = 80
    b_q = 65
    m = 15
    
    # Number of pages available for the outer relation block
    block_size = m - 1
    
    # --- Plan 1: P as the outer relation ---
    # Number of passes over the inner relation Q
    num_passes_p_outer = math.ceil(b_p / block_size)
    cost_p_outer = b_p + (num_passes_p_outer * b_q)
    
    # --- Plan 2: Q as the outer relation ---
    # Number of passes over the inner relation P
    num_passes_q_outer = math.ceil(b_q / block_size)
    cost_q_outer = b_q + (num_passes_q_outer * b_p)

    # Determine the minimum cost and print the result
    print("Block Nested Loop Join (BNLJ) Cost Calculation")
    print(f"Pages in P, B(P) = {b_p}")
    print(f"Pages in Q, B(Q) = {b_q}")
    print(f"Memory pages, M = {m}")
    print("-" * 40)
    
    print("Plan 1: P as outer, Q as inner")
    print(f"Cost = B(P) + (ceil(B(P) / (M - 1)) * B(Q))")
    print(f"Cost = {b_p} + (ceil({b_p} / ({m} - 1)) * {b_q})")
    print(f"Cost = {b_p} + ({num_passes_p_outer} * {b_q}) = {int(cost_p_outer)}")
    print("-" * 40)
    
    print("Plan 2: Q as outer, P as inner")
    print(f"Cost = B(Q) + (ceil(B(Q) / (M - 1)) * B(P))")
    print(f"Cost = {b_q} + (ceil({b_q} / ({m} - 1)) * {b_p})")
    print(f"Cost = {b_q} + ({num_passes_q_outer} * {b_p}) = {int(cost_q_outer)}")
    print("-" * 40)
    
    if cost_q_outer < cost_p_outer:
        min_cost = int(cost_q_outer)
        print("Minimum cost is achieved with Q as the outer relation.")
        print(f"\nFinal Equation: {b_q} + {num_passes_q_outer} * {b_p} = {min_cost}")
        print(f"<<<{min_cost}>>>")
    else:
        min_cost = int(cost_p_outer)
        print("Minimum cost is achieved with P as the outer relation.")
        print(f"\nFinal Equation: {b_p} + {num_passes_p_outer} * {b_q} = {min_cost}")
        print(f"<<<{min_cost}>>>")

calculate_bnlj_cost()