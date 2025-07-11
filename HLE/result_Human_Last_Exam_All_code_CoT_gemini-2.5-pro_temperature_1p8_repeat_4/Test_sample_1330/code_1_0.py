import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # --- Given Parameters ---
    # Number of pages in relation P
    B_P = 80
    # Number of pages in relation Q
    B_Q = 65
    # Number of available memory buffer pages
    M = 15

    print("--- Problem Setup ---")
    print(f"Pages in Relation P, B(P) = {B_P}")
    print(f"Pages in Relation Q, B(Q) = {B_Q}")
    print(f"Available Memory Pages, M = {M}\n")

    print("--- BNLJ Cost Minimization ---")
    print("The Block Nested Loop Join cost formula is: B(outer) + ceil(B(outer) / (M - 2)) * B(inner)")
    print("To minimize I/O cost, the smaller relation is chosen as the outer relation.\n")

    # Determine which relation is smaller to set as the outer relation
    if B_P < B_Q:
        B_outer = B_P
        B_inner = B_Q
        outer_relation_name = "P"
    else:
        B_outer = B_Q
        B_inner = B_P
        outer_relation_name = "Q"
    
    print(f"Comparing relation sizes: B(P)={B_P}, B(Q)={B_Q}.")
    print(f"Since B({outer_relation_name}) is smaller, it will be the outer relation.\n")

    # --- Cost Calculation ---
    print("--- Step-by-Step Calculation ---")

    # Step 1: I/O cost to read the outer relation
    cost_read_outer = B_outer
    print(f"1. Cost to read the outer relation ({outer_relation_name}) once: {cost_read_outer} I/Os")

    # Step 2: Calculate the number of chunks for the outer relation
    # We use M-2 pages for the outer relation block (1 for inner, 1 for output)
    buffer_for_outer = M - 2
    num_chunks = math.ceil(B_outer / buffer_for_outer)
    print(f"2. Number of outer relation chunks = ceil(B(outer) / (M - 2)) = ceil({B_outer} / ({M} - 2)) = {num_chunks}")

    # Step 3: Calculate the I/O cost to scan the inner relation for each chunk
    cost_scan_inner = num_chunks * B_inner
    print(f"3. Cost to scan the inner relation for each chunk = {num_chunks} * {B_inner} = {cost_scan_inner} I/Os")

    # Step 4: Calculate the total cost
    total_cost = cost_read_outer + cost_scan_inner
    print("\n--- Final Cost Equation ---")
    print(f"Total Cost = B(outer) + (ceil(B(outer) / (M - 2)) * B(inner))")
    print(f"Total Cost = {cost_read_outer} + ({num_chunks} * {B_inner})")
    print(f"Total Cost = {cost_read_outer} + {cost_scan_inner}")
    print(f"Minimum I/O Cost = {total_cost}")


if __name__ == "__main__":
    calculate_bnlj_cost()