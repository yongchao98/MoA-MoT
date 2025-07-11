def solve_history_mcq():
    """
    This function identifies the correct statements regarding Wellington's career
    and prints the corresponding numbers in sorted order.
    """
    # Based on historical analysis, the following statements are correct:
    # 1: Correctly links his Indian commissariat system to improvements in the Peninsular campaign.
    # 6: Correctly identifies his method of integrating local forces as a model for later colonial expansion.
    # 8: Correctly notes the adaptation of "flying columns" from India to the Peninsula and later colonial wars like in Burma.
    correct_options = [1, 6, 8]

    # Sort the numbers in ascending order
    correct_options.sort()

    # Format the output as a comma-separated string
    # The numbers in the final output are 1, 6, and 8.
    final_answer = ",".join(map(str, correct_options))
    
    print(final_answer)

solve_history_mcq()