import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join (BNLJ).
    """
    # Given parameters from the problem description
    B_P = 80  # Number of pages in relation P
    B_Q = 65  # Number of pages in relation Q
    M = 15    # Number of available memory buffer pages

    # --- Step 1: Determine the optimal outer and inner relations ---
    # To minimize the BNLJ cost, the smaller relation should be chosen as the outer relation.
    # Cost = B(outer) + ceil(B(outer) / (M - 1)) * B(inner)
    if B_P < B_Q:
        outer_relation_name = "P"
        inner_relation_name = "Q"
        B_outer = B_P
        B_inner = B_Q
    else:
        outer_relation_name = "Q"
        inner_relation_name = "P"
        B_outer = B_Q
        B_inner = B_P
    
    # --- Step 2: Calculate the cost with the optimal configuration ---
    # Number of available buffer pages for the outer relation's block
    buffer_for_outer = M - 1
    
    # Number of times we need to loop over the inner relation
    num_outer_loops = math.ceil(B_outer / buffer_for_outer)
    
    # Total I/O cost
    total_cost = B_outer + (num_outer_loops * B_inner)
    
    # --- Step 3: Print the detailed calculation ---
    print("### Minimum BNLJ Cost Calculation ###")
    print(f"Pages in relation P, B(P) = {B_P}")
    print(f"Pages in relation Q, B(Q) = {B_Q}")
    print(f"Memory buffer pages, M = {M}")
    print("\n--- Strategy ---")
    print("To minimize the join cost, the smaller relation is selected as the outer relation.")
    print(f"Since B({outer_relation_name}) = {B_outer} is less than B({inner_relation_name}) = {B_inner}, we choose {outer_relation_name} as the outer relation.")

    print("\n--- Cost Calculation Equation ---")
    print(f"Cost = B({outer_relation_name}) + (ceil(B({outer_relation_name}) / (M - 1))) * B({inner_relation_name})")
    
    print("\n--- Step-by-Step Execution ---")
    # Output each number in the final equation
    print(f"1. Cost to read the outer relation ({outer_relation_name}): {B_outer} I/Os")
    print(f"2. Number of loops over the inner relation = ceil({B_outer} / ({M} - 1)) = ceil({B_outer} / {buffer_for_outer}) = {num_outer_loops}")
    print(f"3. Cost to scan inner relation ({inner_relation_name}) for all loops = {num_outer_loops} * {B_inner} = {num_outer_loops * B_inner} I/Os")
    print("\n--- Final Cost ---")
    print(f"Total Cost = {B_outer} (read outer) + {num_outer_loops * B_inner} (read inner)")
    print(f"Total Minimum Cost = {total_cost} I/O operations")

# Execute the function
calculate_bnlj_cost()
<<<465>>>