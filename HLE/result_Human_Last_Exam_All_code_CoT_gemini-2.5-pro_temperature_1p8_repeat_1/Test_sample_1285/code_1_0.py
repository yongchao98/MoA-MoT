def solve_problem():
    """
    This script calculates the answer for part (b) of the problem.
    The maximum number of roots of R_t in ]0, 1[ is given by the formula nCr(t, 2),
    which is t * (t - 1) / 2.
    """
    
    # Part (a): State the formula
    # The derivation shows the maximum number of roots is given by the binomial coefficient "t choose 2".
    formula_str_a = "t * (t - 1) / 2"
    print("(a) The expression for the maximum number of roots is: " + formula_str_a)

    # Part (b): Calculate for t = 5
    t = 5
    
    # Calculate the result
    numerator_val_1 = t
    numerator_val_2 = t - 1
    denominator = 2
    max_roots = (numerator_val_1 * numerator_val_2) // denominator
    
    # Print the equation and the result
    print(f"(b) For t = {t}, the calculation is:")
    print(f"{numerator_val_1} * ({t} - 1) / {denominator} = {max_roots}")
    
    # This is the final answer in the format required by the problem prompt.
    # The final answer will be enclosed in <<<>>> at the end of the entire response.

solve_problem()
