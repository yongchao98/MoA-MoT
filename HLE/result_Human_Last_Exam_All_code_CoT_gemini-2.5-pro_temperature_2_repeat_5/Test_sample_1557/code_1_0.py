def solve_wellington_mcq():
    """
    This function identifies the correct statements regarding the Duke of Wellington's career
    and prints the corresponding numbers in sorted order, separated by a comma.

    The correct statements identified through historical analysis are:
    1: Linking the Indian commissariat system to improvements in the Peninsular Campaign.
    6: Stating his methods of using local auxiliary forces became standard colonial practice.
    8: Confirming the adaptation of 'flying columns' from India to the Peninsula and Burma.
    """
    
    # List of numbers for the correct statements
    correct_options = [1, 6, 8]
    
    # Sort the list to ensure the output is in numerical order
    correct_options.sort()
    
    # Convert numbers to strings and join them with a comma
    final_answer = ",".join(map(str, correct_options))
    
    # Print the final result
    print("The correct options are listed below:")
    print(final_answer)

solve_wellington_mcq()