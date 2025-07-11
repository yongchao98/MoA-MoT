def solve_wellington_mcq():
    """
    This function identifies and prints the correct options for the given historical question.
    """
    # Based on historical analysis, the correct statements are 1, 6, and 8.
    correct_options = [1, 6, 8]

    # The options are already in sorted order.
    # The final output should be the numbers separated by a comma.
    answer_string = ",".join(map(str, correct_options))

    print("The correct options listed in sorted order are:")
    print(answer_string)

solve_wellington_mcq()