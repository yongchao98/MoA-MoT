def solve():
    """
    Analyzes the provided mathematical claim and proof, identifies incorrect lines,
    and determines if the claim is correct.
    """
    
    # List of line numbers with incorrect statements
    incorrect_lines = [3, 7]
    
    # The claim is determined to be wrong
    claim_is_correct = "Wrong"
    
    # Format the output string as required
    # The instruction "output each number in the final equation" is interpreted
    # as ensuring the numbers of the incorrect lines are printed.
    line_numbers_str = ", ".join(map(str, incorrect_lines))
    
    print(f"[{line_numbers_str}] {claim_is_correct}")

solve()