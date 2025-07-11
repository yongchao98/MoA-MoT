import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join (BNLJ).
    """
    # Given parameters
    B_P = 80  # Pages in relation P
    B_Q = 65  # Pages in relation Q
    M = 15    # Memory buffer pages

    # In BNLJ, M-2 pages are available for the outer relation's block
    buffer_for_outer = M - 2

    print("Step 1: Define the initial parameters.")
    print(f"Pages in relation P, B(P) = {B_P}")
    print(f"Pages in relation Q, B(Q) = {B_Q}")
    print(f"Available memory pages, M = {M}")
    print(f"Buffer pages for one block of the outer relation = M - 2 = {M} - 2 = {buffer_for_outer}\n")

    # --- Scenario 1: P is the outer relation ---
    print("Step 2: Calculate the cost with P as the outer relation.")
    # Number of blocks for relation P
    num_blocks_p = math.ceil(B_P / buffer_for_outer)
    # Total cost for this scenario
    cost_p_outer = B_P + (num_blocks_p * B_Q)
    
    print("The formula is: Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
    print(f"Number of blocks of P = ceil({B_P} / {buffer_for_outer}) = {num_blocks_p}")
    print(f"Cost = {B_P} + ({num_blocks_p} * {B_Q})")
    cost_inner_scans = num_blocks_p * B_Q
    print(f"     = {B_P} + {cost_inner_scans}")
    print(f"     = {cost_p_outer} I/O operations\n")


    # --- Scenario 2: Q is the outer relation ---
    print("Step 3: Calculate the cost with Q as the outer relation.")
    # Number of blocks for relation Q
    num_blocks_q = math.ceil(B_Q / buffer_for_outer)
    # Total cost for this scenario
    cost_q_outer = B_Q + (num_blocks_q * B_P)

    print("The formula is: Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
    print(f"Number of blocks of Q = ceil({B_Q} / {buffer_for_outer}) = {num_blocks_q}")
    print(f"Cost = {B_Q} + ({num_blocks_q} * {B_P})")
    cost_inner_scans_2 = num_blocks_q * B_P
    print(f"     = {B_Q} + {cost_inner_scans_2}")
    print(f"     = {cost_q_outer} I/O operations\n")


    # --- Conclusion ---
    print("Step 4: Determine the minimum cost.")
    min_cost = min(cost_p_outer, cost_q_outer)
    print(f"Comparing the two costs, the minimum is {min_cost} I/O operations.")
    print("This is achieved when Q, the smaller relation, is used as the outer relation.")

if __name__ == "__main__":
    calculate_bnlj_cost()
    final_answer = 465
    print(f"<<<{final_answer}>>>")