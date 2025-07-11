def solve_history_mcq():
    """
    This function identifies the correct statements about the Duke of Wellington's career.
    
    The correct statements are:
    1: Wellington's mobile commissariat system in India was adapted for the Peninsular Campaign with great success.
    6: His method of integrating local auxiliary forces in India became a standard colonial practice.
    8: His use of 'flying columns' was adapted from India for the Peninsular War and later colonial campaigns like Burma.
    """
    
    # List of the numbers of the correct options
    correct_options = [1, 6, 8]
    
    # Sort the list in ascending order
    correct_options.sort()
    
    # Format the output as a string with numbers separated by commas
    # The final print statement is constructed to clearly show the equation of how we get the final answer.
    print(f"{correct_options[0]},{correct_options[1]},{correct_options[2]}")

solve_history_mcq()