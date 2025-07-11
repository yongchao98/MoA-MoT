def solve_rado_problem():
    """
    This function determines the answers to the theoretical questions based on the logical derivation above.
    """

    # Part (a): Rad_2(S-1) = 1. The answer is "Yes".
    answer_a = "Yes"

    # Part (b): Rad_2(2S-2) can be 2. The answer is "yes".
    # The value is 1 for S=1 and 2 for S>1. An expression for this is floor(S/2) + 1.
    answer_b_yes_no = "yes"
    answer_b_expr = "floor(S/2) + 1"

    # Part (c): Rad_2(2S-1) for even S is 3.
    answer_c = "3"

    # Format the final answer string as requested.
    # The instruction to "output each number in the final equation" is interpreted as
    # providing the numerical answers and the derived expression, which contains numbers.
    final_answer = f"(a) {answer_a}; (b) {answer_b_yes_no} {answer_b_expr}; (c) {answer_c}"

    # Print the final answer in the required wrapper format.
    print(f"<<<{final_answer}>>>")

solve_rado_problem()