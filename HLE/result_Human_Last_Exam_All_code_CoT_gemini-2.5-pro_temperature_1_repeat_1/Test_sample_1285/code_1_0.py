import math

def solve_wronskian_roots():
    """
    This function determines the maximum number of roots for the Wronskian
    determinant Rt as described in the problem. It provides the general formula
    and calculates the specific value for t=5.
    """
    
    # Part (a): The general expression for the maximum number of roots.
    # From mathematical theory, this is given by the binomial coefficient "t choose 2".
    answer_a_expr = "t*(t-1)/2"

    # Part (b): The specific case where t = 5.
    t = 5
    
    # Calculate the maximum number of roots for t=5 using the formula.
    # math.comb(n, k) calculates "n choose k".
    max_roots = math.comb(t, 2)
    
    # Per the instructions, the output for part (b) should show the numbers
    # used in the final calculation.
    answer_b_calc = f"({t} * ({t} - 1)) / 2 = {max_roots}"
    
    # Print the final answer in the specified format: (a) [expression]; (b) [expression].
    print(f"(a) {answer_a_expr}; (b) {answer_b_calc}")

solve_wronskian_roots()