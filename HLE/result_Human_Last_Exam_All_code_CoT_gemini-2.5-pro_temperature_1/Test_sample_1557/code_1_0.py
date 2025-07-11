def solve_history_mcq():
    """
    This function identifies and prints the correct options from the provided list.
    The correct options are determined to be 1, 6, and 8 based on historical analysis.
    """
    correct_options = [1, 6, 8]
    # Sort the options in ascending order
    correct_options.sort()
    # Format the output as a comma-separated string
    result = ",".join(map(str, correct_options))
    print(result)

solve_history_mcq()