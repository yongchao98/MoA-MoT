import math

def calculate_bnlj_cost():
    """
    Calculates the minimum cost of a Block Nested Loop Join.
    """
    # Step 1: Define the given parameters.
    B_P = 80  # Number of pages for relation P
    B_Q = 65  # Number of pages for relation Q
    M = 15    # Number of available memory buffer pages

    # Step 2: Determine the optimal outer and inner relations.
    # To minimize cost, the smaller relation is chosen as the outer relation.
    if B_P < B_Q:
        B_outer = B_P
        B_inner = B_Q
        outer_relation_name = "P"
        inner_relation_name = "Q"
    else:
        B_outer = B_Q
        B_inner = B_P
        outer_relation_name = "Q"
        inner_relation_name = "P"

    # Step 3: Calculate the BNLJ cost.
    # The cost formula is: B(outer) + ceil(B(outer) / (M - 2)) * B(inner)
    
    # Buffer pages for the outer relation blocks
    buffer_for_outer = M - 2
    
    # Number of times the inner relation must be scanned
    num_inner_scans = math.ceil(B_outer / buffer_for_outer)
    
    # Total I/O cost
    total_cost = B_outer + (num_inner_scans * B_inner)

    # --- Output the calculation step-by-step ---
    print("Calculating the minimum Block Nested Loop Join (BNLJ) cost:")
    print("-" * 60)
    print(f"Given: B(P) = {B_P}, B(Q) = {B_Q}, Memory Pages (M) = {M}")
    print(f"To minimize cost, the smaller relation ({outer_relation_name}) is chosen as the outer relation.")
    print("-" * 60)
    print("BNLJ Cost Formula: B(outer) + ceil(B(outer) / (M - 2)) * B(inner)")
    print("\nSubstituting the values:")
    print(f"Cost = B({outer_relation_name}) + ceil(B({outer_relation_name}) / ({M} - 2)) * B({inner_relation_name})")
    print(f"Cost = {B_outer} + ceil({B_outer} / {buffer_for_outer}) * {B_inner}")
    print(f"Cost = {B_outer} + {num_inner_scans} * {B_inner}")
    print(f"Cost = {B_outer} + {num_inner_scans * B_inner}")
    print(f"Minimum BNLJ Cost = {total_cost}")


if __name__ == "__main__":
    calculate_bnlj_cost()
