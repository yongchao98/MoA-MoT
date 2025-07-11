def solve_history_mcq():
    """
    This function evaluates the historical statements about the Duke of Wellington and identifies the correct ones.

    The correct statements are:
    1: Wellington's mobile commissariat system developed in India was successfully adapted for the Peninsular Campaign.
    6: Wellington's methods of integrating local auxiliary forces, first tested in India, became standard practice in British colonial expansion.
    8: Wellington's system of deploying flying columns, developed in India, was adapted for the Peninsular War and later colonial campaigns.

    The function will print the numbers of these correct statements in sorted order, separated by a comma.
    """
    correct_options = [1, 6, 8]
    
    # Sort the numbers just in case they were not added in order
    correct_options.sort()
    
    # Format the output string
    result = ",".join(map(str, correct_options))
    
    print(result)

solve_history_mcq()