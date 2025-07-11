import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join (BNLJ).
    """
    # Given parameters
    B_P = 80  # Number of pages in relation P
    B_Q = 65  # Number of pages in relation Q
    M = 15    # Number of available memory buffer pages

    # The number of buffer pages available for holding blocks of the outer relation.
    # We reserve one page for reading blocks of the inner relation.
    buffer_for_outer = M - 1

    # --- Scenario 1: P as outer relation, Q as inner relation ---
    num_outer_scans_1 = math.ceil(B_P / buffer_for_outer)
    cost_P_outer = B_P + (num_outer_scans_1 * B_Q)

    # --- Scenario 2: Q as outer relation, P as inner relation ---
    num_outer_scans_2 = math.ceil(B_Q / buffer_for_outer)
    cost_Q_outer = B_Q + (num_outer_scans_2 * B_P)

    # --- Determine the minimum cost and print the optimal scenario ---
    print("Calculating the minimum I/O cost for the Block Nested Loop Join.")
    print("-" * 60)
    print(f"Given Parameters:")
    print(f"  - Pages in P, B(P): {B_P}")
    print(f"  - Pages in Q, B(Q): {B_Q}")
    print(f"  - Memory Buffer Pages, M: {M}")
    print("-" * 60)

    # To minimize the cost, we should choose the smaller relation as the outer relation.
    if cost_P_outer < cost_Q_outer:
        print("Optimal strategy: Use P as the outer relation and Q as the inner relation.\n")
        print("Cost Formula: B(P) + ceil(B(P) / (M - 1)) * B(Q)")
        print(f"Calculation:")
        print(f"Cost = {B_P} + ceil({B_P} / ({M} - 1)) * {B_Q}")
        print(f"Cost = {B_P} + ceil({B_P} / {buffer_for_outer}) * {B_Q}")
        print(f"Cost = {B_P} + {num_outer_scans_1} * {B_Q}")
        print(f"Cost = {B_P} + {num_outer_scans_1 * B_Q}")
        print(f"Cost = {cost_P_outer}")
        min_cost = cost_P_outer
    else:
        print("Optimal strategy: Use Q as the outer relation and P as the inner relation.\n")
        print("Cost Formula: B(Q) + ceil(B(Q) / (M - 1)) * B(P)")
        print(f"Calculation:")
        print(f"Cost = {B_Q} + ceil({B_Q} / ({M} - 1)) * {B_P}")
        print(f"Cost = {B_Q} + ceil({B_Q} / {buffer_for_outer}) * {B_P}")
        # To show the division result before ceiling
        division_result = B_Q / buffer_for_outer
        print(f"Cost = {B_Q} + ceil({division_result:.2f}) * {B_P}")
        print(f"Cost = {B_Q} + {num_outer_scans_2} * {B_P}")
        print(f"Cost = {B_Q} + {num_outer_scans_2 * B_P}")
        print(f"Cost = {cost_Q_outer}")
        min_cost = cost_Q_outer
    
    print("-" * 60)
    print(f"Minimum I/O Cost: {min_cost}")

if __name__ == "__main__":
    calculate_bnlj_cost()
    print("<<<465>>>")
