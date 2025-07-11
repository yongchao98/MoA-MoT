def solve_string_theory_problem():
    """
    Solves the three-part problem about sl(2) strings in diagonal harmonics.
    """
    
    # Part a: Find the bi-degree of the terminal polynomial.
    # A string starter P (highest weight vector) with bi-degree (a, b) generates a string
    # of length k+1, where k = a - b. The terminal polynomial has bi-degree (a-k, b+k).
    # Given (a, b) = (4, 3).
    a1, b1 = 4, 3
    k1 = a1 - b1
    # The expression for the terminal bi-degree is (a1 - k1, b1 + k1).
    final_a = a1 - k1
    final_b = b1 + k1
    
    # We format the string to show the calculation as requested.
    answer_a_expr = f"({a1} - {k1}, {b1} + {k1})"
    
    # Part b: Find the condition for a polynomial of bi-degree (a, b) to be a starter.
    # For a polynomial to be a starter (highest weight vector) of a full sl(2) string,
    # its weight, a-b, must be a non-negative integer. So, a >= b.
    answer_b_expr = "a >= b"
    
    # Part c: Check if a polynomial of bi-degree (5, 2) can be a starter.
    # We apply the condition from part b. Let a=5, b=2.
    # Since 5 >= 2, the primary condition is satisfied.
    # The existence of such a polynomial in DH_n for some n (e.g., n>=7) is known,
    # so it is possible.
    a3, b3 = 5, 2
    if a3 >= b3:
        answer_c_expr = "Yes"
    else:
        answer_c_expr = "No"
        
    # Print the final answer in the required format "a) [Expression] b) [Expression] c) [Yes/No]".
    # Note: For part a, we show the expression that evaluates to the final answer.
    print(f"a) {answer_a_expr} b) {answer_b_expr} c) {answer_c_expr}")

solve_string_theory_problem()