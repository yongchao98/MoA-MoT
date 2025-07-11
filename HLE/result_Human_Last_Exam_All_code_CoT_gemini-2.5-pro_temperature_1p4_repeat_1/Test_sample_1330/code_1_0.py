import math

def calculate_bnlj_cost():
    """
    Calculates the minimum cost for a Block Nested Loop Join (BNLJ) given two relations P and Q.
    """
    # Given parameters
    B_P = 80  # Pages in relation P
    B_Q = 65  # Pages in relation Q
    M = 15    # Available memory buffer pages

    # For BNLJ, we use M-2 pages for the outer relation block, reserving one page for
    # the inner relation's input and one for the output buffer.
    buffer_for_outer = M - 2

    # --- Scenario 1: P is the outer relation, Q is the inner ---
    num_blocks_P = math.ceil(B_P / buffer_for_outer)
    cost_P_outer = B_P + (num_blocks_P * B_Q)

    # --- Scenario 2: Q is the outer relation, P is the inner ---
    num_blocks_Q = math.ceil(B_Q / buffer_for_outer)
    cost_Q_outer = B_Q + (num_blocks_Q * B_P)

    # --- Determine the minimum cost and print the result ---
    print("Block Nested Loop Join (BNLJ) Cost Calculation")
    print("------------------------------------------------")
    print(f"Given Relations: P ({B_P} pages), Q ({B_Q} pages)")
    print(f"Memory Buffer: {M} pages")
    print(f"Buffer available for outer relation block: M - 2 = {buffer_for_outer} pages")
    print("\nThere are two possible execution plans:\n")

    print(f"1. P as the outer relation:")
    print(f"   Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
    print(f"   Cost = {B_P} + (ceil({B_P} / {buffer_for_outer}) * {B_Q})")
    print(f"   Cost = {B_P} + ({num_blocks_P} * {B_Q}) = {B_P} + {num_blocks_P * B_Q} = {cost_P_outer} I/Os\n")

    print(f"2. Q as the outer relation:")
    print(f"   Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
    print(f"   Cost = {B_Q} + (ceil({B_Q} / {buffer_for_outer}) * {B_P})")
    print(f"   Cost = {B_Q} + ({num_blocks_Q} * {B_P}) = {B_Q} + {num_blocks_Q * B_P} = {cost_Q_outer} I/Os\n")
    
    # Conclusion
    if cost_Q_outer < cost_P_outer:
        print("Conclusion: The minimum cost is achieved when Q is the outer relation.")
        print("Final Minimum Cost Calculation:")
        print(f"Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
        print(f"     = {B_Q} + (ceil({B_Q} / ({M} - 2)) * {B_P})")
        print(f"     = {B_Q} + (ceil({B_Q} / {buffer_for_outer}) * {B_P})")
        print(f"     = {B_Q} + ({num_blocks_Q} * {B_P})")
        print(f"     = {B_Q} + {num_blocks_Q * B_P}")
        print(f"     = {cost_Q_outer}")
        final_answer = cost_Q_outer
    else:
        print("Conclusion: The minimum cost is achieved when P is the outer relation.")
        print("Final Minimum Cost Calculation:")
        print(f"Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
        print(f"     = {B_P} + (ceil({B_P} / ({M} - 2)) * {B_Q})")
        print(f"     = {B_P} + (ceil({B_P} / {buffer_for_outer}) * {B_Q})")
        print(f"     = {B_P} + ({num_blocks_P} * {B_Q})")
        print(f"     = {B_P} + {num_blocks_P * B_Q}")
        print(f"     = {cost_P_outer}")
        final_answer = cost_P_outer

    # This part is for the final answer submission as per the instructions
    # <<<answer>>>
    return final_answer

# Execute the function and get the final answer
min_cost_value = calculate_bnlj_cost()
print(f"<<<{min_cost_value}>>>")
