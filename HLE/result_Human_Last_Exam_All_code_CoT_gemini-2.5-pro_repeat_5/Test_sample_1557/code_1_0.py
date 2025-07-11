def solve_wellington_mcq():
    """
    This function identifies and prints the correct options for the given historical question.
    Based on a historical analysis of the Duke of Wellington's career, the following statements
    were identified as accurate:

    1: Wellington's mobile commissariat system from India was adapted for the Peninsular Campaign.
    6: Wellington's method of integrating local auxiliary forces in India became a standard colonial practice.
    8: Wellington's tactic of using flying columns, developed in India, was adapted for the Peninsular War and other colonial campaigns.

    The function will print the numbers of these correct options in sorted order, separated by a comma.
    """
    correct_options = [1, 6, 8]
    
    # Sort the list to ensure the output is in ascending order
    correct_options.sort()
    
    # Format the output as a comma-separated string
    output_string = ",".join(map(str, correct_options))
    
    print(output_string)

solve_wellington_mcq()