def solve_rado_questions():
    """
    This function analyzes the theoretical questions about Rado numbers and prints the formatted answer.
    """
    # Analysis results
    # (a) The statement is true. Rad_2(S-1) = 1.
    answer_a = "Yes"

    # (b) Rad_2(2S-2) can equal 2 (this occurs when S > 1).
    answer_b_can = "yes"
    answer_b_expression = 2

    # (c) For c = 2S-1 and even S, the Rado number is 2S-1.
    # The prompt requests to output each number in the final equation.
    # The expression for the value is "2*S - 1". The numbers are 2 and 1.
    answer_c_num1 = 2
    answer_c_num2 = 1

    # Print the answers in the specified format.
    print(f"(a) {answer_a}; (b) {answer_b_can} [{answer_b_expression}]; (c) {answer_c_num1}*S - {answer_c_num2}")

solve_rado_questions()