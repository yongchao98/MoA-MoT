def solve_trivia():
    """
    This function stores and prints the answer to the historical trivia question.
    """
    # The feature that was removed.
    feature = "A tramway"
    
    # The year mentioned in the question.
    year = 1950

    # The prompt asks to "output each number in the final equation",
    # but there is no equation. We will just print the relevant year from the prompt.
    print(f"The unique architectural feature removed from the Piazza della Rotonda around the year {year} was:")
    print(feature)

solve_trivia()