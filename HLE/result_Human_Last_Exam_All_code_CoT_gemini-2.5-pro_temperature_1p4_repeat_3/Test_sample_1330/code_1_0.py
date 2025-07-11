import math

def calculate_bnlj_cost():
    """
    Calculates the minimum cost for a Block Nested Loop Join (BNLJ)
    given the page sizes of two relations and the available memory buffer size.
    """
    # Given parameters
    B_P = 80  # Number of pages for relation P
    B_Q = 65  # Number of pages for relation Q
    M = 15    # Number of available memory buffer pages

    print("--- Calculating Minimum Block Nested Loop Join Cost ---")
    print(f"Relation P pages B(P) = {B_P}")
    print(f"Relation Q pages B(Q) = {B_Q}")
    print(f"Memory buffer pages M = {M}")
    print("-" * 50)

    # --- Scenario 1: P is the outer relation ---
    print("Scenario 1: P as outer relation, Q as inner relation")
    # Number of times we need to read the entire inner relation
    num_outer_blocks_p = math.ceil(B_P / (M - 2))
    # Total I/O cost for this scenario
    cost_p_outer = B_P + num_outer_blocks_p * B_Q

    print(f"Cost = B(P) + (ceil(B(P) / (M - 2))) * B(Q)")
    print(f"Cost = {B_P} + (ceil({B_P} / ({M} - 2))) * {B_Q}")
    print(f"Cost = {B_P} + {num_outer_blocks_p} * {B_Q}")
    print(f"Cost = {B_P} + {num_outer_blocks_p * B_Q}")
    print(f"Total cost with P as outer: {cost_p_outer} I/Os")
    print("-" * 50)

    # --- Scenario 2: Q is the outer relation ---
    print("Scenario 2: Q as outer relation, P as inner relation")
    # Number of times we need to read the entire inner relation
    num_outer_blocks_q = math.ceil(B_Q / (M - 2))
    # Total I/O cost for this scenario
    cost_q_outer = B_Q + num_outer_blocks_q * B_P

    print(f"Cost = B(Q) + (ceil(B(Q) / (M - 2))) * B(P)")
    print(f"Cost = {B_Q} + (ceil({B_Q} / ({M} - 2))) * {B_P}")
    print(f"Cost = {B_Q} + {num_outer_blocks_q} * {B_P}")
    print(f"Cost = {B_Q} + {num_outer_blocks_q * B_P}")
    print(f"Total cost with Q as outer: {cost_q_outer} I/Os")
    print("-" * 50)

    # --- Conclusion ---
    min_cost = min(cost_p_outer, cost_q_outer)
    if min_cost == cost_q_outer:
        print("The minimum cost is achieved when Q is the outer relation.")
    else:
        print("The minimum cost is achieved when P is the outer relation.")
    
    print(f"\nThe minimum cost for the join is {min_cost} I/O operations.")

if __name__ == '__main__':
    calculate_bnlj_cost()
    print("<<<465>>>")
