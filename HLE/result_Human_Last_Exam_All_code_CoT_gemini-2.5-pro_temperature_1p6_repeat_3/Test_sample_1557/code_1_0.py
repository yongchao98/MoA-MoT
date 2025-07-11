def solve_history_mcq():
    """
    Analyzes the historical statements about the Duke of Wellington and identifies the correct ones.
    The correct statements are determined to be 1, 6, and 8 based on historical analysis.
    This function will print the numbers of the correct options, sorted and separated by a comma.
    """
    correct_options = [1, 6, 8]
    # Sort the numbers just in case they were not identified in order
    correct_options.sort()
    # Format the output as a comma-separated string
    print(",".join(map(str, correct_options)))

solve_history_mcq()