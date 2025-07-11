import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Relation P has 80 pages
    bp = 80
    # Relation Q has 65 pages
    bq = 65
    # The system has 15 buffer pages
    m = 15

    # To minimize cost, the smaller relation should be the outer relation.
    if bq < bp:
        b_outer = bq
        b_inner = bp
        outer_relation_name = "Q"
        inner_relation_name = "P"
    else:
        b_outer = bp
        b_inner = bq
        outer_relation_name = "P"
        inner_relation_name = "Q"

    # Number of buffer pages available for the block of the outer relation
    # is M-2 (1 for the inner relation page, 1 for the output buffer).
    buffer_for_outer = m - 2

    # Number of blocks the outer relation is divided into.
    # We must scan the entire inner relation for each of these blocks.
    num_blocks = math.ceil(b_outer / buffer_for_outer)

    # Total I/O cost = cost to read outer relation + (cost to read inner relation for each outer block)
    total_cost = b_outer + (num_blocks * b_inner)

    print("--- Block Nested Loop Join (BNLJ) Cost Calculation ---")
    print(f"Pages in P, B(P) = {bp}")
    print(f"Pages in Q, B(Q) = {bq}")
    print(f"Buffer pages, M = {m}\n")
    print(f"To minimize cost, we choose the smaller relation ({outer_relation_name}) as the outer relation.")

    print("\n--- Cost Calculation Steps ---")
    print(f"Cost = B({outer_relation_name}) + (ceil(B({outer_relation_name}) / (M - 2)) * B({inner_relation_name}))")
    print(f"Cost = {b_outer} + (ceil({b_outer} / ({m} - 2)) * {b_inner})")
    print(f"Cost = {b_outer} + (ceil({b_outer} / {buffer_for_outer}) * {b_inner})")
    print(f"Cost = {b_outer} + ({int(num_blocks)} * {b_inner})")
    print(f"Cost = {b_outer} + {int(num_blocks * b_inner)}")
    print(f"Final Minimum Cost = {int(total_cost)} I/O operations")

if __name__ == "__main__":
    calculate_bnlj_cost()