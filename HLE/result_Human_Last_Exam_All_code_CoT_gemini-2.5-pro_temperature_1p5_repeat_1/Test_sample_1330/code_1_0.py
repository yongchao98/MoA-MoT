import math

def calculate_bnlj_cost():
    """
    Calculates the minimum cost for a Block Nested Loop Join (BNLJ)
    and prints the step-by-step calculation.
    """
    # Given parameters from the problem description
    B_P = 80  # Number of pages in relation P
    B_Q = 65  # Number of pages in relation Q
    M = 15    # Number of available memory buffer pages

    # To minimize BNLJ cost, the smaller relation should be the outer relation.
    # Here, B(Q) < B(P), so Q is the outer relation and P is the inner relation.
    B_outer = B_Q
    B_inner = B_P
    outer_relation_name = "Q"
    inner_relation_name = "P"

    # The block nested loop join algorithm uses M-2 pages for the outer relation's block,
    # 1 page for the inner relation's tuples, and 1 for the output buffer.
    buffer_for_outer = M - 2

    # Calculate the number of blocks for the outer relation.
    # This is the number of times we will need to scan the entire inner relation.
    num_blocks_outer = math.ceil(B_outer / buffer_for_outer)

    # The total cost is the cost of reading the outer relation once, plus
    # the cost of reading the inner relation for each block of the outer relation.
    total_cost = B_outer + (num_blocks_outer * B_inner)

    print("To minimize the Block Nested Loop Join (BNLJ) cost, the smaller relation is chosen as the outer relation.")
    print(f"Pages in P: {B_P}, Pages in Q: {B_Q}. Since {B_Q} < {B_P}, Q is the outer relation.\n")

    print("The cost is calculated using the formula:")
    print("Cost = B(Outer) + (ceil(B(Outer) / (M - 2)) * B(Inner))\n")

    print("Substituting the values into the formula:")
    print(f"Cost = B({outer_relation_name}) + (ceil(B({outer_relation_name}) / (M - 2)) * B({inner_relation_name}))")
    print(f"Cost = {B_outer} + (ceil({B_outer} / ({M} - 2)) * {B_inner})")
    print(f"Cost = {B_outer} + (ceil({B_outer} / {buffer_for_outer}) * {B_inner})")
    print(f"Cost = {B_outer} + ({num_blocks_outer} * {B_inner})")
    print(f"Cost = {B_outer} + {num_blocks_outer * B_inner}")
    print(f"Minimum Cost = {total_cost} I/O operations")

calculate_bnlj_cost()
<<<465>>>