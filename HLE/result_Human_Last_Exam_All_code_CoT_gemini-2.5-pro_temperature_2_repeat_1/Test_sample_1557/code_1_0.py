def find_correct_statements():
    """
    This function analyzes the historical statements about the Duke of Wellington
    and identifies the correct ones based on established historical facts.

    The correct statements are:
    1: Wellington's logistical innovations in India were foundational for the Peninsular Campaign.
    6: His method of integrating local auxiliary forces in India became standard colonial practice.
    8: His tactic of using flying columns, developed in India, was adapted for other wars.
    """

    # The numbers of the historically correct statements
    correct_options = [1, 6, 8]

    # The list is already sorted, but we ensure it for good practice
    correct_options.sort()

    # The final output should be the numbers of the correct options,
    # sorted and separated by a comma.
    # The format requirement is "output each number in the final equation"
    # which we interpret as printing the final answer clearly.
    final_answer_string = ",".join(map(str, correct_options))
    
    print(f"<<<{final_answer_string}>>>")

find_correct_statements()