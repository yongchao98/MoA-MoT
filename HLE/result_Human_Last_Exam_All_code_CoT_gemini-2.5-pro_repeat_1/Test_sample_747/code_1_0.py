def solve_task():
    """
    Analyzes the provided mathematical proof and identifies the incorrect lines and the validity of the claim.
    """
    # List of line numbers with incorrect statements
    incorrect_lines = [3, 7, 8]

    # Verdict on the claim
    claim_is_correct = "Wrong"

    # Format the final answer string as per the user's request
    # The final equation/answer is the string "[3,7,8] Wrong"
    # We output each component of this string as requested.
    # (Although a simple print("[3,7,8] Wrong") would suffice,
    # the prompt mentioned outputting each number in the final equation,
    # so we construct the string programmatically.)
    
    result_str = f"[{','.join(map(str, incorrect_lines))}] {claim_is_correct}"
    print(result_str)

solve_task()