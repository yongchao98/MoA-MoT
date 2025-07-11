def solve_wronskian_roots():
    """
    Calculates the maximum number of roots for the Wronskian determinant R_t.
    
    (a) Finds the general expression for the maximum number of roots.
    (b) Computes this maximum number for t = 5.
    """
    
    # Part (a): The general expression for the maximum number of roots in ]0, 1[
    # is given by the formula C(t, 2) = t*(t-1)/2.
    expression_a = "t*(t-1)/2"

    # Part (b): Calculate the maximum number of roots for t = 5.
    # We apply the formula derived in part (a).
    t = 5
    
    # The numbers in the final equation are t, t-1, and 2.
    numerator_term_1 = t
    numerator_term_2 = t - 1
    denominator = 2
    
    # The equation to calculate the result is (t * (t-1)) / 2
    result_b = (numerator_term_1 * numerator_term_2) // denominator

    # The final output is printed in the specified format.
    print(f"(a) {expression_a}; (b) {result_b}.")

# Execute the function to print the solution.
solve_wronskian_roots()