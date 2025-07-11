import math

def calculate_bnlj_cost():
    """
    Calculates the minimum cost of a Block Nested Loop Join (BNLJ).
    """
    # Step 1: Identify Given Parameters
    B_P = 80  # Number of pages for relation P
    B_Q = 65  # Number of pages for relation Q
    M = 15    # Number of available memory buffer pages

    # Step 2: Understand BNLJ memory usage
    # We use M-2 pages for the outer relation's block.
    M_usable = M - 2

    print("--- Calculating Minimum BNLJ Cost ---")
    print(f"Pages in Relation P, B(P): {B_P}")
    print(f"Pages in Relation Q, B(Q): {B_Q}")
    print(f"Available Memory Pages, M: {M}")
    print(f"Usable Buffer for Outer Relation, M-2: {M_usable}\n")
    print("The cost formula is: B(outer) + ceil(B(outer) / (M - 2)) * B(inner)\n")

    # Step 3 & 4: Calculate costs for both scenarios

    # --- Scenario 1: P as Outer Relation ---
    print("Scenario 1: P is Outer, Q is Inner")
    # Number of blocks P is divided into
    num_blocks_P = math.ceil(B_P / M_usable)
    cost_P_outer = B_P + num_blocks_P * B_Q
    
    print(f"Number of blocks for P = ceil({B_P} / {M_usable}) = {num_blocks_P}")
    print(f"Cost = B(P) + (Blocks of P) * B(Q)")
    print(f"Equation: {B_P} + {num_blocks_P} * {B_Q} = {B_P} + {num_blocks_P * B_Q} = {cost_P_outer} I/Os\n")

    # --- Scenario 2: Q as Outer Relation ---
    print("Scenario 2: Q is Outer, P is Inner")
    # Number of blocks Q is divided into
    num_blocks_Q = math.ceil(B_Q / M_usable)
    cost_Q_outer = B_Q + num_blocks_Q * B_P
    
    print(f"Number of blocks for Q = ceil({B_Q} / {M_usable}) = {num_blocks_Q}")
    print(f"Cost = B(Q) + (Blocks of Q) * B(P)")
    print(f"Equation: {B_Q} + {num_blocks_Q} * {B_P} = {B_Q} + {num_blocks_Q * B_P} = {cost_Q_outer} I/Os\n")
    
    # Step 5: Determine the minimum cost
    min_cost = min(cost_P_outer, cost_Q_outer)
    
    print("--- Conclusion ---")
    print("To minimize cost, we choose the relation that requires fewer blocks as the outer relation.")
    if min_cost == cost_Q_outer:
        print("Choosing Q as the outer relation is cheaper.")
    else:
        print("Choosing P as the outer relation is cheaper.")
        
    print(f"\nThe minimum cost for the Block Nested Loop Join is {min_cost} I/O operations.")
    return min_cost

if __name__ == "__main__":
    calculate_bnlj_cost()