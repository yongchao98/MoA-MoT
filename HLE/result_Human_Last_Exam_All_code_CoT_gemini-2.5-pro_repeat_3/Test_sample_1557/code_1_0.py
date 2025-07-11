def solve_wellington_mcq():
    """
    This function identifies and prints the correct options for the given historical question.
    """
    # Based on the historical analysis, the following statements are correct:
    # 1: Wellington's commissariat system from India was adapted for the Peninsular Campaign.
    # 6: Methods of integrating local forces from India became standard colonial practice.
    # 8: The use of flying columns from India was adapted for the Peninsula and Burma.
    correct_options = [1, 6, 8]

    # The options are already in sorted order.
    # We will format the output as a comma-separated string.
    # The prompt requests to "output each number in the final equation!", which we interpret
    # as clearly printing each number of the final answer set.
    
    answer_string = ",".join(map(str, correct_options))
    print(f"<<<{answer_string}>>>")

solve_wellington_mcq()