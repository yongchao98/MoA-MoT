import math

def calculate_bnlj_cost():
    """
    Calculates the minimum I/O cost for a Block Nested Loop Join.
    """
    # Step 1: Define the given parameters
    b_p = 80  # Pages in relation P
    b_q = 65  # Pages in relation Q
    m = 15    # Available memory buffer pages

    # In BNLJ, M-2 pages are available for the outer relation's block.
    # One page is for scanning the inner relation, one is for the output buffer.
    buffer_for_outer = m - 2

    # Step 2: Calculate cost for both scenarios to find the minimum.
    # Generally, making the smaller relation the outer one is optimal.

    # Scenario 1: P is outer, Q is inner
    num_loops_p_outer = math.ceil(b_p / buffer_for_outer)
    cost_p_outer = b_p + num_loops_p_outer * b_q

    # Scenario 2: Q is outer, P is inner
    num_loops_q_outer = math.ceil(b_q / buffer_for_outer)
    cost_q_outer = b_q + num_loops_q_outer * b_p

    # Step 3: Determine the minimum cost and print the detailed calculation for that case.
    print("--- Block Nested Loop Join (BNLJ) Cost Calculation ---")
    print(f"Pages in P, B(P) = {b_p}")
    print(f"Pages in Q, B(Q) = {b_q}")
    print(f"Memory Buffer Pages, M = {m}")
    print(f"Buffer for Outer Relation Block = M - 2 = {m} - 2 = {buffer_for_outer}")
    print("\nTo minimize I/O cost, we choose the smaller relation as the outer relation for the join.")
    print(f"Since B(Q) < B(P) ({b_q} < {b_p}), we will use Q as the outer relation.")

    # Select the optimal case (Q as outer)
    min_cost = cost_q_outer
    b_outer = b_q
    b_inner = b_p
    num_loops = num_loops_q_outer

    print("\n--- Minimum Cost Calculation (Q as Outer) ---")
    print(f"Cost Formula: B(Outer) + ceil(B(Outer) / (M - 2)) * B(Inner)")
    # Print the equation with initial values
    print(f"Cost = {b_outer} + ceil({b_outer} / {buffer_for_outer}) * {b_inner}")
    # Print the equation with the ceiling calculation result
    print(f"Cost = {b_outer} + {num_loops} * {b_inner}")
    # Print the equation after the multiplication
    print(f"Cost = {b_outer} + {num_loops * b_inner}")
    # Print the final result
    print(f"Cost = {min_cost}")

    print(f"\nThe minimum cost to execute the query is {min_cost} I/O operations.")


if __name__ == '__main__':
    calculate_bnlj_cost()
