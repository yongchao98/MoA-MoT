def solve():
    """
    Analyzes the provided mathematical claim and proof, identifies the incorrect lines
    in the proof, and determines if the claim is correct.
    """
    # List of line numbers with incorrect statements
    incorrect_lines = [3, 7, 8]

    # Evaluation of the claim's correctness
    claim_is_correct = "Wrong"

    # The problem asks to print the list of incorrect lines and the evaluation of the claim.
    # The format is [line1, line2, ...] Right/Wrong
    print(f"[{', '.join(map(str, incorrect_lines))}] {claim_is_correct}")

solve()