def solve_wellington_mcq():
    """
    This function identifies and prints the correct options for the multiple-choice question
    regarding the Duke of Wellington's career.
    """
    # After analyzing each statement, the correct options are identified.
    # 1: Correct. Indian commissariat experience was adapted for the Peninsula.
    # 3: Correct. His divisional structure became a standard for the British Army.
    # 6: Correct. Integrating local auxiliary forces became standard colonial practice.
    # 8: Correct. The use of 'flying columns' was adapted from India for other conflicts.
    correct_options = [1, 3, 6, 8]
    
    # Sorting the options in ascending order as requested.
    correct_options.sort()
    
    # Formatting the output as a comma-separated string.
    answer = ",".join(map(str, correct_options))
    
    print("The correct options listed in sorted order are:")
    print(f"<<<{answer}>>>")

solve_wellington_mcq()