def solve():
    """
    Calculates the maximum number of roots for R_t and provides the solution in the specified format.
    """
    t = 5

    # (a) The expression for the maximum number of roots in terms of t.
    expression_a = "t(t-1)/2"

    # (b) The maximum number for t = 5.
    t_val = 5
    result_b = t_val * (t_val - 1) // 2
    
    # Format the calculation for part (b) as requested.
    calc_b = f"{t_val}({t_val}-1)/2 = {result_b}"

    # Print the final answer in the format (a) [expression]; (b) [expression].
    print(f"(a) {expression_a}; (b) {calc_b}")

solve()