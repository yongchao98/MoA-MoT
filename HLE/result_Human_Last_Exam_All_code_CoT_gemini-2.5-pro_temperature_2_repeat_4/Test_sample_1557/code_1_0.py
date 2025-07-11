def solve_task():
    """
    This function analyzes the historical statements about the Duke of Wellington and prints the numbers of the correct statements.
    
    The correct statements identified through historical analysis are:
    1: Regarding the adaptation of the Indian mobile commissariat for the Peninsular Campaign.
    6: Regarding the integration of local auxiliary forces becoming standard practice.
    8: Regarding the adaptation of flying columns for use in the Peninsular War and later in Burma.
    
    The numbers of the correct statements, sorted, are 1, 6, and 8.
    """
    
    # The correct options determined by historical analysis
    correct_options = [1, 6, 8]
    
    # Sorting the numbers in ascending order
    correct_options.sort()
    
    # Formatting the output as a string with comma separators
    # The prompt asks to "output each number in the final equation"
    # This interprets as printing the final list of numbers clearly.
    final_answer = ",".join(map(str, correct_options))
    
    print(final_answer)

solve_task()