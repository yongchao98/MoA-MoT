def solve_history_mcq():
    """
    This function identifies and prints the correct options for the given MCQ.
    The correct options, based on historical analysis, are 1, 6, and 8.
    """
    correct_options = [1, 6, 8]
    
    # The options are already in sorted order.
    # The output format requires each number from the final answer to be printed, separated by a comma.
    
    # Convert numbers to strings to join them
    options_as_strings = []
    for option in correct_options:
        options_as_strings.append(str(option))
        
    final_answer = ",".join(options_as_strings)
    
    print(final_answer)

solve_history_mcq()