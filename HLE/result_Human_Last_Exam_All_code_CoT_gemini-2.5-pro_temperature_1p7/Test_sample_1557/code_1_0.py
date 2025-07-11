def solve_wellington_mcq():
    """
    This function identifies the correct statements regarding Wellington's career
    and prints the corresponding numbers in sorted order.
    """
    # Based on historical analysis, the following statements are correct:
    # 1: Wellington's mobile commissariat in India was adapted for the Peninsular Campaign.
    # 6: Methods of integrating local auxiliary forces became standard colonial practice.
    # 8: The use of 'flying columns' was adapted for the Peninsula and later colonial wars.
    correct_options = [1, 6, 8]

    # Sort the options to ensure they are in numerical order.
    correct_options.sort()

    # Format the output as a comma-separated string.
    # The final equation requested is the list of correct numbers.
    # The problem asks to "output each number in the final equation",
    # so we will print them separated by commas.
    final_answer = ",".join(map(str, correct_options))
    
    print(final_answer)

solve_wellington_mcq()