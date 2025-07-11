import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    B_P = 80
    B_Q = 65
    M = 15

    # Case 1: P is the outer relation
    cost_P_outer = B_P + math.ceil(B_P / (M - 2)) * B_Q

    # Case 2: Q is the outer relation
    cost_Q_outer = B_Q + math.ceil(B_Q / (M - 2)) * B_P

    # Determine the minimum cost and print the calculation
    if cost_P_outer < cost_Q_outer:
        print("Minimum cost is achieved with P as the outer relation.")
        num_chunks = math.ceil(B_P / (M - 2))
        final_cost = cost_P_outer
        
        print(f"Cost = B(P) + ceil(B(P) / (M - 2)) * B(Q)")
        print(f"Cost = {B_P} + ceil({B_P} / ({M} - 2)) * {B_Q}")
        print(f"Cost = {B_P} + ceil({B_P} / {M - 2}) * {B_Q}")
        print(f"Cost = {B_P} + {num_chunks} * {B_Q}")
        print(f"Cost = {B_P} + {num_chunks * B_Q}")
        print(f"Cost = {final_cost}")

    else:
        print("Minimum cost is achieved with Q as the outer relation.")
        num_chunks = math.ceil(B_Q / (M - 2))
        final_cost = cost_Q_outer
        
        print(f"Cost = B(Q) + ceil(B(Q) / (M - 2)) * B(P)")
        print(f"Cost = {B_Q} + ceil({B_Q} / ({M} - 2)) * {B_P}")
        print(f"Cost = {B_Q} + ceil({B_Q} / {M - 2}) * {B_P}")
        print(f"Cost = {B_Q} + {num_chunks} * {B_P}")
        print(f"Cost = {B_Q} + {num_chunks * B_P}")
        print(f"Cost = {final_cost}")

calculate_bnlj_cost()