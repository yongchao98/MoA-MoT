import math

def calculate_bnlj_cost():
    """
    Calculates the minimum cost for a Block Nested Loop Join.
    """
    # Given parameters
    b_p = 80  # Pages in relation P
    b_q = 65  # Pages in relation Q
    m = 15    # Available memory buffer pages

    # In a Block Nested Loop Join, one page is for the inner relation's input,
    # and one is for the output buffer. The rest are for the outer relation's block.
    buffer_for_outer = m - 2

    # To minimize the cost, we must choose the smaller relation as the outer relation.
    # Here, B(Q) < B(P), so Q will be the outer relation and P the inner.
    b_outer = b_q
    b_inner = b_p
    outer_relation_name = "Q"
    inner_relation_name = "P"

    # Step 1: Calculate the number of blocks the outer relation (Q) is split into.
    # This is the number of times we will have to scan the entire inner relation.
    num_blocks_outer = math.ceil(b_outer / buffer_for_outer)

    # Step 2: Calculate the total I/O cost.
    # Cost = (Pages read from outer) + (Number of outer blocks * Pages in inner)
    total_cost = b_outer + (num_blocks_outer * b_inner)
    
    # --- Output the explanation and calculation ---
    print(f"To minimize the join cost, the smaller relation is chosen as the outer relation.")
    print(f"Given B(P) = {b_p} and B(Q) = {b_q}, relation Q is smaller and is chosen as the outer relation.")
    print("")
    print("The BNLJ cost is calculated as: Cost = B(outer) + (ceil(B(outer) / (M - 2)) * B(inner))")
    print("")
    
    # Show the calculation for the number of blocks
    print("1. Calculate the number of blocks for the outer relation (Q):")
    print(f"   num_blocks_Q = ceil(B(Q) / (M - 2))")
    print(f"   num_blocks_Q = ceil({b_outer} / ({m} - 2))")
    print(f"   num_blocks_Q = ceil({b_outer} / {buffer_for_outer}) = {num_blocks_outer}")
    print("")

    # Show the final cost calculation with all the numbers
    print("2. Calculate the total I/O cost:")
    print(f"   Cost = B(Q) + (num_blocks_Q * B(P))")
    print(f"   Cost = {b_outer} + ({num_blocks_outer} * {b_inner})")
    cost_part_2 = num_blocks_outer * b_inner
    print(f"   Cost = {b_outer} + {cost_part_2}")
    print(f"   Cost = {total_cost}")

    # Final answer in the required format
    print(f"\n<<<465>>>")

calculate_bnlj_cost()