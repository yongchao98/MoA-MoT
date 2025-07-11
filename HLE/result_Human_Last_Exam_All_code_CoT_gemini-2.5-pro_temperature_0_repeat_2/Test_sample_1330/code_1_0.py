import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Given parameters
    B_P = 80  # Number of pages for relation P
    B_Q = 65  # Number of pages for relation Q
    M = 15    # Number of available memory buffer pages

    # The BNLJ cost formula is: B(outer) + ceil(B(outer) / (M - 2)) * B(inner)
    # We must choose the smaller relation as the outer relation to minimize cost.
    # Since B(Q) < B(P), we choose Q as the outer relation.

    B_outer = B_Q
    B_inner = B_P
    outer_relation_name = "Q"
    inner_relation_name = "P"

    # Number of pages available for the outer relation's block
    buffer_for_outer = M - 2

    # Number of blocks the outer relation is divided into. This is also the
    # number of times the inner relation must be scanned.
    num_outer_blocks = math.ceil(B_outer / buffer_for_outer)

    # Cost to scan the outer relation once
    cost_read_outer = B_outer

    # Cost to scan the inner relation for each block of the outer relation
    cost_scan_inner = num_outer_blocks * B_inner

    # Total minimum cost
    min_cost = cost_read_outer + cost_scan_inner

    print("Calculating the minimum cost for the Block Nested Loop Join (BNLJ).")
    print("The optimal strategy is to use the smaller relation (Q) as the outer relation.")
    print("\n--- Final Equation ---")
    print(f"Cost = B({outer_relation_name}) + ceil(B({outer_relation_name}) / (M - 2)) * B({inner_relation_name})")
    print(f"Cost = {B_outer} + ceil({B_outer} / ({M} - 2)) * {B_inner}")
    print(f"Cost = {B_outer} + ceil({B_outer} / {buffer_for_outer}) * {B_inner}")
    print(f"Cost = {B_outer} + {num_outer_blocks} * {B_inner}")
    print(f"Cost = {cost_read_outer} + {cost_scan_inner}")
    print(f"Cost = {min_cost}")

calculate_bnlj_cost()