def find_next_term():
    """
    Analyzes the sequence to find the next term.
    """
    # The given sequence
    s = [2, 11, 23, 51, 119]
    
    print("Analyzing the relationship between consecutive terms:")
    print("The pattern appears to be next_term = 2 * current_term + residual.\n")
    
    residuals = []
    # Loop to find the residual for each step
    for i in range(len(s) - 1):
        current_term = s[i]
        next_term = s[i+1]
        residual = next_term - 2 * current_term
        residuals.append(residual)
        print(f"{next_term} = 2 * {current_term} + {residual}")

    print("\nThe sequence of residuals is:", residuals)
    
    # We observe that the pattern in the residuals starts from the second residual.
    # Let's analyze the sub-sequence of residuals starting from the second one.
    sub_residuals = residuals[1:]
    print("\nAnalyzing the pattern in the residuals, starting from the second one:", sub_residuals)
    
    # Calculate the differences between these residuals
    residual_diffs = []
    for i in range(len(sub_residuals) - 1):
        diff = sub_residuals[i+1] - sub_residuals[i]
        residual_diffs.append(diff)
    
    print(f"The differences between these residuals are {residual_diffs[0]} and {residual_diffs[1]}.")
    print("The sequence of differences [4, 12] is a geometric progression with a ratio of 3 (12 / 4 = 3).")
    
    # Calculate the next difference and the next residual
    next_diff = residual_diffs[-1] * 3
    print(f"\nThe next difference should be 12 * 3 = {next_diff}.")
    
    next_residual = sub_residuals[-1] + next_diff
    print(f"The next residual should be 17 + {next_diff} = {next_residual}.")
    
    # Calculate the final term in the main sequence
    last_known_term = s[-1]
    next_term_in_sequence = 2 * last_known_term + next_residual
    
    print("\nTherefore, the next term in the sequence is calculated as follows:")
    print(f"{int(next_term_in_sequence)} = 2 * {last_known_term} + {next_residual}")

find_next_term()