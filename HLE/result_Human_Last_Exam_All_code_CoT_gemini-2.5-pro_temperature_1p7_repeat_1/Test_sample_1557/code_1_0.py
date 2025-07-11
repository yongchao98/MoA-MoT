def find_correct_statements():
    """
    This function identifies the correct statements regarding the Duke of Wellington's career.
    Based on historical analysis, the following options are correct:

    1: Wellington's mobile commissariat system developed during the Siege of Seringapatam (1799) was later successfully adapted for the Peninsular Campaign.
    6: Wellington's methods of integrating local auxiliary forces, first tested in India, became standard practice in British colonial expansion during the 1830s and 1840s.
    8: Wellington's system of deploying flying columns, developed in India, was successfully adapted for use in both the Peninsular War and later colonial campaigns in Burma.

    The function will return these numbers as a comma-separated string.
    """
    
    # List of correct option numbers
    correct_options = [1, 6, 8]
    
    # Sort the list to ensure the output is in numerical order
    correct_options.sort()
    
    # Format the numbers into the final string as requested.
    # The statement "output each number in the final equation" is interpreted
    # as clearly printing each of the correct option numbers.
    final_answer = ", ".join(map(str, correct_options))
    
    print(final_answer)

find_correct_statements()