def solve():
    """
    Solves the problem by providing the expression for part (a) and calculating the result for part (b).
    """

    # Part (a): Provide the general expression for the maximum number of roots.
    # The sum symbol is represented as 'sum', and min as 'min'.
    part_a_expression = "(sum_{i=1 to t} (k_i + l_i)) - t*(min(k_1,...,k_t) + min(l_1,...,l_t)) - t*(t-1)/2"

    # Part (b): Calculate the value for t=5.
    t = 5
    # As explained in the plan, under the assumption of a canonical choice of exponents
    # (the Bernstein basis), the expression simplifies to t*(t-1)/2.
    
    # We show the calculation step-by-step.
    numerator = t * (t - 1)
    result = numerator / 2
    
    # Format the output as requested.
    part_b_expression = f"{t} * ({t}-1) / 2 = {int(result)}"

    print(f"(a) {part_a_expression}; (b) {part_b_expression}.")

solve()