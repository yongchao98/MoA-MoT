def solve_set_problem():
    """
    Calculates the number of sets (S1, S2, S3, S4) satisfying the given conditions.
    """
    # Number of choices for placing each element into the regions defined by the set chain.
    choices_for_1 = 1  # 1 must be in S1
    choices_for_2 = 2  # 2 must be in S2 (i.e., in S1 or S2\S1)
    choices_for_3 = 3  # 3 must be in S3 (i.e., in S1, S2\S1, or S3\S2)
    choices_for_4 = 5  # 4 can be anywhere
    choices_for_5 = 5  # 5 can be anywhere

    # The total number of ways is the product of the choices for each element.
    total_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

    # Print the explanation and the calculation
    print("The number of choices for each element is determined by the constraints:")
    print(f"Choices for element 1: {choices_for_1}")
    print(f"Choices for element 2: {choices_for_2}")
    print(f"Choices for element 3: {choices_for_3}")
    print(f"Choices for element 4: {choices_for_4}")
    print(f"Choices for element 5: {choices_for_5}")
    print("\nThe total number of sets is the product of these choices:")
    print(f"Total = {choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_sets}")

solve_set_problem()