def solve_set_problem():
    """
    This function calculates the number of set tuples (S1, S2, S3, S4) satisfying the given conditions.
    """
    print("Let U = {1, 2, 3, 4, 5} be the universal set.")
    print("We are looking for the number of tuples of sets (S1, S2, S3, S4) such that:")
    print("1. S1 is a subset of S2, which is a subset of S3, which is a subset of S4, which is a subset of U.")
    print("2. i is in Si for i = 1, 2, 3.")
    print("\nThis problem can be solved by considering the possible placements for each element of U.\n")
    print("The chain of subsets S1 <= S2 <= S3 <= S4 <= U defines 5 disjoint regions for the elements:")
    print("  - Region 1: in S1")
    print("  - Region 2: in S2 but not S1")
    print("  - Region 3: in S3 but not S2")
    print("  - Region 4: in S4 but not S3")
    print("  - Region 5: in U but not S4")
    print("\nWe count the number of choices for each element based on the constraints:\n")

    # Number of choices for element 1
    # Constraint: 1 in S1. This means 1 must be in Region 1.
    choices_for_1 = 1
    print(f"Choices for element 1 (must be in S1): {choices_for_1}")

    # Number of choices for element 2
    # Constraint: 2 in S2. This means 2 can be in Region 1 or Region 2.
    choices_for_2 = 2
    print(f"Choices for element 2 (must be in S2): {choices_for_2}")

    # Number of choices for element 3
    # Constraint: 3 in S3. This means 3 can be in Region 1, 2, or 3.
    choices_for_3 = 3
    print(f"Choices for element 3 (must be in S3): {choices_for_3}")

    # Number of choices for element 4
    # No constraint. Can be in any of the 5 regions.
    choices_for_4 = 5
    print(f"Choices for element 4 (no constraint): {choices_for_4}")

    # Number of choices for element 5
    # No constraint. Can be in any of the 5 regions.
    choices_for_5 = 5
    print(f"Choices for element 5 (no constraint): {choices_for_5}")

    # The total number of valid tuples is the product of these choices.
    total_count = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5
    
    print("\nThe total number of sets is the product of the choices for each element.")
    print(f"The calculation is: {choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_count}")
    
solve_set_problem()