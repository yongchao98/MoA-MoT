def solve():
    """
    This function determines the answers to the three parts of the question
    and prints them in the specified format.
    """
    answer_a = "True"
    answer_b = "True"
    
    # The expression for part (c) is binom(s-1, 2) + (s-1)*(t-1).
    # This can be written in standard Python arithmetic as:
    expression_c = "(s-1)*(s-2)//2 + (s-1)*(t-1)"
    
    # Format the final output string as requested.
    # The instruction "output each number in the final equation" is interpreted as
    # presenting the formula with its constituent numbers (1, 2, etc.) visible.
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"
    
    print(final_answer)

solve()