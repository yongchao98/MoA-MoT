def solve_set_problem():
    """
    Calculates the number of sets (S1, S2, S3, S4) satisfying
    S1 subset S2 subset S3 subset S4 subset {1,2,3,4,5} and i in Si for i=1,2,3.
    """
    # Let U = {1, 2, 3, 4, 5}. The chain of subsets partitions U into 5 disjoint regions.
    # We count the number of placement choices for each element of U.

    # Choices for element 1: 1 must be in S1.
    choices_for_1 = 1

    # Choices for element 2: 2 must be in S2. Can be in S1 or S2\S1.
    choices_for_2 = 2

    # Choices for element 3: 3 must be in S3. Can be in S1, S2\S1, or S3\S2.
    choices_for_3 = 3

    # Choices for element 4: No restrictions. Can be in any of the 5 regions.
    choices_for_4 = 5

    # Choices for element 5: No restrictions. Can be in any of the 5 regions.
    choices_for_5 = 5

    # The total number of ways is the product of the choices for each element.
    total_ways = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

    print("The total number of sets is the product of the number of placement choices for each element.")
    print(f"Number of choices for element 1: {choices_for_1}")
    print(f"Number of choices for element 2: {choices_for_2}")
    print(f"Number of choices for element 3: {choices_for_3}")
    print(f"Number of choices for element 4: {choices_for_4}")
    print(f"Number of choices for element 5: {choices_for_5}")
    print("The final calculation is:")
    print(f"{choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_ways}")

solve_set_problem()