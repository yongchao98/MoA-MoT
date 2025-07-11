def solve():
    """
    Analyzes the mathematical claim and proof to identify incorrect lines and evaluate the claim.
    """
    # List of line numbers with incorrect statements
    incorrect_lines = [3, 7, 8]

    # Verdict on the claim itself
    claim_is_correct = "Wrong"

    # The final answer format is a list of numbers and a word.
    # The problem asks to print the numbers in the final equation,
    # which we interpret as printing the elements of our list.
    print(f"[{', '.join(map(str, incorrect_lines))}] {claim_is_correct}")

solve()