def solve_history_mcq():
    """
    This function identifies the correct statements about the Duke of Wellington's career
    and prints them in the required format.
    """

    # Based on historical analysis, the following statements are correct:
    # 1: Correctly links his Indian commissariat experience to the Peninsular War.
    # 3: Correctly identifies the influence of his organizational structures on the post-1815 British Army.
    # 6: Correctly notes that his method of integrating local forces became standard colonial practice.
    # 8: Correctly states that the use of "flying columns" was a tactic developed in India and adapted elsewhere.
    correct_options = [1, 3, 6, 8]

    # The options are already in sorted order.
    # We will now format the output as a string of numbers separated by a comma.
    answer = ",".join(map(str, correct_options))

    # Print the final answer in the specified format.
    print(f"<<<{answer}>>>")

solve_history_mcq()